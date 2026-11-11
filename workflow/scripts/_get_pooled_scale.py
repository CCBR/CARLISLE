#!/usr/bin/env python3
"""
Compute the scaling factor for a pooled control bedgraph by summing read counts
across matching control replicates from alignment_stats.tsv.

Column and scale logic:
  LIBRARY  -> lib_factor / sum(nreads_genome) across matching replicates
               nreads_genome is dedup_nreads_genome (dupstatus=dedup) or
               no_dedup_nreads_genome (dupstatus=no_dedup).
               lib_factor = tiered power-of-10 constant derived from median(nreads_genome) across ALL
               samples — mirrors _make_library_norm_table.R exactly, so pooled-control
               bedgraphs land on the same normalisation scale as individual replicates.
  SPIKEIN  -> spikein_scale / sum(no_dedup_nreads_spikein) across matching replicates
               (fragment-length + mapq filtered, duplicates intentionally NOT removed;
                matches the spike-in counts used in the individual-replicate bam2bg rule)
  NONE     -> scale = 1

Output: a single floating-point scaling factor printed to stdout.
"""

import argparse
import csv
import re
import statistics
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--align_stats",
        required=True,
        help="Path to alignment_stats.tsv",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--control_sample",
        required=False,
        help=(
            "Literal control sample name (without replicate suffix). "
            "Replicates are matched as '^<escaped_control_sample>_[0-9]+$'."
        ),
    )
    group.add_argument(
        "--sample_pattern",
        required=False,
        help="Regex pattern to match control replicate sample_name values "
        "(e.g. '^HN6_IgG_rabbit_negative_control_[0-9]+$')",
    )
    parser.add_argument(
        "--norm_method",
        required=True,
        choices=["LIBRARY", "SPIKEIN", "NONE"],
        help="Normalization method",
    )
    parser.add_argument(
        "--dupstatus",
        required=False,
        default="dedup",
        choices=["dedup", "no_dedup"],
        help=(
            "Duplication status of the pooled bedgraph being scaled. "
            "Used only for LIBRARY normalization to select the correct read-count column "
            "(dedup_nreads_genome vs no_dedup_nreads_genome). Default: dedup"
        ),
    )
    parser.add_argument(
        "--spikein_scale",
        type=float,
        default=1.0,
        help="Scaling constant (numerator) for SPIKEIN mode. Ignored for LIBRARY mode. Default: 1.0",
    )
    return parser.parse_args()


def _lib_factor(median_reads: float) -> float:
    """Return the library size factor, mirroring the tiered if/else logic in
    _make_library_norm_table.R.

    Uses strict > comparisons identical to R (e.g. median == 1000 falls into
    the 1e2 bucket, not 1e3). Clamps to 1e1 for medians <= 100. Caps at 1e8
    for very large libraries, matching the R script's highest tier.
    """
    if median_reads > 1e8:
        return 1e8
    elif median_reads > 1e7:
        return 1e7
    elif median_reads > 1e6:
        return 1e6
    elif median_reads > 1e5:
        return 1e5
    elif median_reads > 1e4:
        return 1e4
    elif median_reads > 1e3:
        return 1e3
    elif median_reads > 1e2:
        return 1e2
    else:
        return 1e1


def main():
    args = parse_args()

    if args.norm_method == "NONE":
        print(1)
        return

    # Choose read-count column
    if args.norm_method == "LIBRARY":
        # dupstatus determines which read-count column to use, matching
        # _make_library_norm_table.R which runs once per dedup_type
        column = (
            "dedup_nreads_genome"
            if args.dupstatus == "dedup"
            else "no_dedup_nreads_genome"
        )
    else:
        # SPIKEIN: always use no_dedup spike-in counts (duplicates not removed)
        column = "no_dedup_nreads_spikein"

    if args.control_sample is not None:
        pattern_text = rf"^{re.escape(args.control_sample)}_[0-9]+$"
    else:
        pattern_text = args.sample_pattern
    pattern = re.compile(pattern_text)
    all_values: list = []  # all samples — used to compute lib_factor median (LIBRARY only)
    ctrl_total = 0.0
    ctrl_count = 0

    with open(args.align_stats, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames  # triggers header read
        if not fieldnames:
            sys.exit(f"ERROR: {args.align_stats} is empty or has no header row")
        if column not in fieldnames:
            sys.exit(
                f"ERROR: column '{column}' not found in {args.align_stats}.\n"
                f"Available columns: {fieldnames}"
            )
        for row in reader:
            try:
                val = float(row[column])
            except (ValueError, KeyError) as e:
                sys.exit(f"ERROR: could not parse column '{column}' in row {row}: {e}")
            if args.norm_method == "LIBRARY":
                all_values.append(val)  # collect every sample for median
            if pattern.search(row["sample_name"]):
                ctrl_total += val
                ctrl_count += 1

    if ctrl_count == 0:
        sys.exit(
            f"ERROR: no rows matched pattern '{pattern_text}' in {args.align_stats}"
        )

    if ctrl_total == 0:
        sys.exit(
            f"ERROR: sum of '{column}' across matched replicates is zero. "
            "Cannot compute scaling factor."
        )

    if args.norm_method == "LIBRARY":
        if not all_values:
            sys.exit(f"ERROR: {args.align_stats} has no data rows")
        factor = _lib_factor(statistics.median(all_values))
        print(factor / ctrl_total)
    else:
        # SPIKEIN
        print(args.spikein_scale / ctrl_total)


if __name__ == "__main__":
    main()
