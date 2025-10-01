#!/usr/bin/env python

# MIT License
# (c) 2025 — adapted from logic by @chris-cheshire
# Upstream reference: https://github.com/nf-core/cutandrun/blob/6e1125d4fee4ea7c8b70ed836bb0e92a89e3305f/bin/find_unique_reads.py

"""
BAM linear-duplicate filtering with pysam

Overview
--------
This utility filters BAM alignments with two possible modes:

1) Default (linear-dedup mode)
   - For *paired-end*: deduplication key is `chrom + strand(read1) + pos(read1)`.
     The mate's position is intentionally ignored (to reflect linear-amplification duplicate logic).
     Only primary, mapped alignments are considered; read1 alone defines keys.
     A fragment-length filter is applied using |tlen| from read1 when available.
     For each key, the read with the highest MAPQ is kept (ties resolved by last seen).
     All primary alignments with kept QNAMEs (both mates) are written to the output BAM.
   - For *single-end*: key is `chrom + strand + pos`, highest MAPQ per key is kept.

2) Filter-only mode (`--no-linear-dedup`)
   - No deduplication is performed.
   - Only primary, mapped reads passing MAPQ are retained.
   - For *paired-end*, reads must also be **properly paired** (`is_proper_pair`) and pass the fragment-length filter (|tlen| ≤ `--fraglen` when available).
   - For *single-end*, only MAPQ and primary mapping filters apply.

Common behavior
---------------
- Secondary and supplementary alignments are excluded from consideration and output.
- Output is BAM with original header.
- A plain-text metrics file summarizes counts. The labels reflect the selected mode.

Assumptions & caveats
---------------------
- Input BAM can be name- or coordinate-sorted; the tool reads sequentially (`until_eof=True`).
- Fragment-length filter uses `abs(template_length)`; if 0 or undefined (e.g., discordant/split situations),
  length filtering is skipped in both modes.
- Memory usage grows with number of unique start sites in dedup mode (hash map of unique keys).
- PCR/UMI-aware workflows are out of scope here (no UMI parsing).

Examples
--------
# Deduplicate + filters (default)
python dedup_bam.py --bam_in in.bam --bam_out out.bam --metrics_path stats.txt --minmapq 30 --fraglen 1000

# Filter-only: primary + proper pair (PE) + MAPQ + fragment length; no dedup
python dedup_bam.py --bam_in in.bam --bam_out out.bam --metrics_path stats.txt --no-linear-dedup --minmapq 30 --fraglen 1000
"""

import argparse
import pysam

Description = (
    "Filter BAM using pysam. Default: linear-amplification dedup (PE key = chrom+strand+pos1), "
    "keeping highest MAPQ per key. With --no-linear-dedup: only filter primary (and proper-pair for PE), "
    "MAPQ, and fragment length."
)

parser = argparse.ArgumentParser(description=Description)
parser.add_argument("--bam_in", required=True, help="Input BAM")
parser.add_argument("--bam_out", required=True, help="Output BAM")
parser.add_argument("--metrics_path", required=True, help="Metrics text output")
parser.add_argument("--fraglen", type=int, default=int(1e6),
                    help="Max allowed fragment length for paired-end (default: 1e6)")
parser.add_argument("--minmapq", type=int, default=20,
                    help="Minimum mapping quality (default: 20)")
parser.add_argument("--no-linear-dedup", action="store_true",
                    help="Disable linear deduplication; apply only primary/proper-pair (PE), MAPQ, and fragment-length filters")
args = parser.parse_args()


def is_primary_mapped(r):
    return (not r.is_unmapped) and (not r.is_secondary) and (not r.is_supplementary)


def strand_char(r):
    return "+" if not r.is_reverse else "-"


def filter_only_mode():
    """
    No deduplication: write reads that are primary (and proper pair for PE),
    pass MAPQ, and pass fragment-length filter for PE.
    """
    i = 0            # reads considered
    kept = 0         # reads written

    with pysam.AlignmentFile(args.bam_in, "rb") as bam_in, \
         pysam.AlignmentFile(args.bam_out, "wb", header=bam_in.header) as bam_out:

        for r in bam_in.fetch(until_eof=True):
            if not is_primary_mapped(r):
                continue
            if r.mapping_quality < args.minmapq:
                continue

            if r.is_paired:
                # Require proper pairing in filter-only mode
                if not r.is_proper_pair:
                    continue
                tlen = abs(r.template_length)
                if tlen != 0 and tlen > args.fraglen:
                    continue

            # Count after passing all filters
            i += 1
            bam_out.write(r)
            kept += 1

    return i, kept


def dedup_mode():
    """
    Linear dedup:
    - PE dedup key: chrom + strand(read1) + pos(read1); ignore mate position.
    - SE dedup key: chrom + strand + pos.
    - Keep highest-MAPQ per key (store QNAMEs), then write all primary reads whose QNAMEs are kept.
    """
    # First pass: decide which QNAMEs to keep
    paired_best = {}  # key -> (qname, mapq)
    single_best = {}

    considered = 0   # number of alignments that entered dedup logic

    with pysam.AlignmentFile(args.bam_in, "rb") as bam_in:
        for r in bam_in.fetch(until_eof=True):
            if not is_primary_mapped(r):
                continue
            if r.mapping_quality < args.minmapq:
                continue

            ref = bam_in.get_reference_name(r.reference_id)
            pos = r.reference_start  # 0-based
            if r.is_paired:
                # Only read1 defines the key; mate position is ignored by design.
                if not r.is_read1:
                    continue

                tlen = abs(r.template_length)
                if tlen != 0 and tlen > args.fraglen:
                    continue

                key = f"{ref}{strand_char(r)}{pos}"
                considered += 1
                if (key not in paired_best) or (paired_best[key][1] < r.mapping_quality):
                    paired_best[key] = (r.query_name, r.mapping_quality)
            else:
                key = f"{ref}{strand_char(r)}{pos}"
                considered += 1
                if (key not in single_best) or (single_best[key][1] < r.mapping_quality):
                    single_best[key] = (r.query_name, r.mapping_quality)

    kept_qnames = set(q for q, _ in paired_best.values()) | set(q for q, _ in single_best.values())
    unique_sites = len(kept_qnames)  # count of unique start-site winners (QNAMEs)

    # Second pass: write all primary reads for the kept QNAMEs
    written_reads = 0
    with pysam.AlignmentFile(args.bam_in, "rb") as bam_in, \
         pysam.AlignmentFile(args.bam_out, "wb", header=bam_in.header) as bam_out:

        for r in bam_in.fetch(until_eof=True):
            if not is_primary_mapped(r):
                continue
            if r.query_name in kept_qnames:
                bam_out.write(r)
                written_reads += 1

    # Metrics for dedup mode:
    # - considered: number of reads that participated in dedup (read1 for PE and all SE)
    # - unique_sites: number of unique alignments retained (winners)
    return considered, unique_sites


def main():
    if args.no_linear_dedup:
        # FILTER-ONLY MODE
        i, kept = filter_only_mode()
        removed = i - kept
        removed_pct = round((removed / i * 100.0), 2) if i else 0.0

        report = []
        report.append("FILTER-ONLY METRICS (no linear deduplication)")
        report.append(f"Mode\tfilter-only")
        report.append(f"Min MAPQ\t{args.minmapq}")
        report.append(f"Max fragment length (PE)\t{args.fraglen}")
        report.append(f"Reads passing filters\t{kept}")
        report.append(f"Reads failing filters (n)\t{removed}")
        report.append(f"Reads failing filters (%)\t{removed_pct}")

        with open(args.metrics_path, "w") as f:
            f.write("\n".join(report))
    else:
        # LINEAR DEDUP MODE
        considered, unique_kept = dedup_mode()
        dup_removed = considered - unique_kept
        dup_pct = round((dup_removed / considered * 100.0), 2) if considered else 0.0

        report = []
        report.append("LINEAR AMPLIFICATION DUPLICATION METRICS")
        report.append(f"Mode\tlinear-dedup")
        report.append(f"Min MAPQ\t{args.minmapq}")
        report.append(f"Max fragment length (PE)\t{args.fraglen}")
        report.append(f"Reads before dedup (considered)\t{considered}")
        report.append(f"LA duplicates removed (n)\t{dup_removed}")
        report.append(f"LA duplicates removed (%)\t{dup_pct}")
        report.append(f"Unique alignments after LA duplicate removal\t{unique_kept}")

        with open(args.metrics_path, "w") as f:
            f.write("\n".join(report))


if __name__ == "__main__":
    main()
