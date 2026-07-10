#!/usr/bin/env python3
"""Merge replicate peak BEDs into treatment-level consensus peak sets.

Consensus policy:
- overlap >= overlap_bp_min counts as support
- keep clusters supported by at least min_replicate_support distinct replicates
- merged interval spans leftest-left/rightest-right
- width guard: if merged width exceeds max_merged_width_bp, re-cluster with strict_overlap_bp
  and if still too wide keep interval but mark width_guard_flag=1
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


@dataclass
class Peak:
    chrom: str
    start: int
    end: int
    replicate: str


def parse_input_spec(spec: str) -> Tuple[str, str]:
    if "::" not in spec:
        raise ValueError(f"Invalid --input value '{spec}'. Expected replicate::path format.")
    replicate, path = spec.split("::", 1)
    if not replicate:
        raise ValueError(f"Invalid --input value '{spec}'. Missing replicate id.")
    if not path:
        raise ValueError(f"Invalid --input value '{spec}'. Missing path.")
    return replicate, path


def load_peaks(input_specs: Sequence[str]) -> Dict[str, List[Peak]]:
    by_chrom: Dict[str, List[Peak]] = defaultdict(list)
    for spec in input_specs:
        replicate, path = parse_input_spec(spec)
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                by_chrom[chrom].append(Peak(chrom=chrom, start=start, end=end, replicate=replicate))
    return by_chrom


def cluster_intervals(peaks: Sequence[Peak], overlap_bp: int) -> List[List[Peak]]:
    if not peaks:
        return []
    sorted_peaks = sorted(peaks, key=lambda p: (p.start, p.end))
    clusters: List[List[Peak]] = []
    current = [sorted_peaks[0]]
    current_end = sorted_peaks[0].end

    for peak in sorted_peaks[1:]:
        # overlap condition with configurable minimum overlap length.
        if peak.start <= current_end - overlap_bp:
            current.append(peak)
            if peak.end > current_end:
                current_end = peak.end
        else:
            clusters.append(current)
            current = [peak]
            current_end = peak.end

    clusters.append(current)
    return clusters


def cluster_to_record(
    cluster: Sequence[Peak],
    min_support: int,
    max_width: int,
    strict_overlap: int,
) -> List[Tuple[str, int, int, int, str, int, int]]:
    reps = sorted({p.replicate for p in cluster})
    if len(reps) < min_support:
        return []

    start = min(p.start for p in cluster)
    end = max(p.end for p in cluster)
    width = end - start

    if width <= max_width:
        return [
            (
                cluster[0].chrom,
                start,
                end,
                len(reps),
                ",".join(reps),
                len(cluster),
                0,
            )
        ]

    # Too wide: re-cluster with stricter overlap and evaluate each subcluster.
    refined_records: List[Tuple[str, int, int, int, str, int, int]] = []
    refined_clusters = cluster_intervals(cluster, strict_overlap)
    for sub in refined_clusters:
        sub_reps = sorted({p.replicate for p in sub})
        if len(sub_reps) < min_support:
            continue
        sub_start = min(p.start for p in sub)
        sub_end = max(p.end for p in sub)
        sub_width = sub_end - sub_start
        refined_records.append(
            (
                sub[0].chrom,
                sub_start,
                sub_end,
                len(sub_reps),
                ",".join(sub_reps),
                len(sub),
                1 if sub_width > max_width else 0,
            )
        )

    return refined_records


def write_records(records: Sequence[Tuple[str, int, int, int, str, int, int]], output_path: str) -> None:
    with open(output_path, "w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(
                "\t".join(
                    [
                        rec[0],
                        str(rec[1]),
                        str(rec[2]),
                        str(rec[3]),
                        rec[4],
                        str(rec[5]),
                        str(rec[6]),
                    ]
                )
                + "\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge replicate peaks into treatment-level consensus sets")
    parser.add_argument("--input", action="append", required=True, help="replicate::path to peak BED")
    parser.add_argument("--output", required=True, help="Output BED path")
    parser.add_argument("--overlap-bp-min", type=int, default=1)
    parser.add_argument("--min-replicate-support", type=int, default=2)
    parser.add_argument("--max-merged-width-bp", type=int, default=10000)
    parser.add_argument("--strict-overlap-bp", type=int, default=50)
    args = parser.parse_args()

    by_chrom = load_peaks(args.input)
    all_records: List[Tuple[str, int, int, int, str, int, int]] = []

    for chrom in sorted(by_chrom):
        clusters = cluster_intervals(by_chrom[chrom], args.overlap_bp_min)
        for cluster in clusters:
            all_records.extend(
                cluster_to_record(
                    cluster,
                    min_support=args.min_replicate_support,
                    max_width=args.max_merged_width_bp,
                    strict_overlap=args.strict_overlap_bp,
                )
            )

    all_records.sort(key=lambda r: (r[0], r[1], r[2]))
    write_records(all_records, args.output)


if __name__ == "__main__":
    main()
