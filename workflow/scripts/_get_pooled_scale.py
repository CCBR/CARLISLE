#!/usr/bin/env python3
"""
Compute the scaling factor for a pooled control bedgraph by summing read counts
across matching control replicates from alignment_stats.tsv.

Column logic:
  LIBRARY  -> sum of dedup_nreads_genome across matching replicates
  SPIKEIN  -> sum of no_dedup_nreads_spikein across matching replicates
               (fragment-length + mapq filtered, duplicates intentionally NOT removed;
                matches the spike-in counts used in the individual-replicate bam2bg rule)
  NONE     -> scale = 1

Output: a single floating-point scaling factor printed to stdout.
"""

import argparse
import csv
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--align_stats",
        required=True,
        help="Path to alignment_stats.tsv",
    )
    parser.add_argument(
        "--sample_pattern",
        required=True,
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
        "--spikein_scale",
        type=float,
        default=1.0,
        help="Scaling constant (numerator). Default: 1.0",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.norm_method == "NONE":
        print(1)
        return

    column = {
        "LIBRARY": "dedup_nreads_genome",
        "SPIKEIN": "no_dedup_nreads_spikein",
    }[args.norm_method]

    pattern = re.compile(args.sample_pattern)
    total = 0.0
    count = 0

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
            if pattern.search(row["sample_name"]):
                try:
                    total += float(row[column])
                    count += 1
                except (ValueError, KeyError) as e:
                    sys.exit(f"ERROR: could not parse value in row {row}: {e}")

    if count == 0:
        sys.exit(
            f"ERROR: no rows matched sample_pattern '{args.sample_pattern}' "
            f"in {args.align_stats}"
        )

    if total == 0:
        sys.exit(
            f"ERROR: sum of '{column}' across matched replicates is zero. "
            "Cannot compute scaling factor."
        )

    print(args.spikein_scale / total)


if __name__ == "__main__":
    main()
