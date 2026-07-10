#!/usr/bin/env python3
"""Aggregate treatment-level HOMER annotation summaries across treatments."""

from __future__ import annotations

import argparse
import os


def parse_input_spec(spec: str):
    if "::" not in spec:
        raise ValueError(f"Invalid --input '{spec}', expected label::path")
    label, path = spec.split("::", 1)
    return label, path


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate treatment-level HOMER summaries")
    parser.add_argument("--input", action="append", required=True, help="label::path")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("source\tannotation\tdistance_to_tss\tnumber_of_peaks\tpercent_of_peaks\ttotal_size_bp\tlog10_p_value\tlog2_ratio\tlogP_enrichment\n")
        for spec in args.input:
            source, path = parse_input_spec(spec)
            if not os.path.exists(path):
                continue
            with open(path, "r", encoding="utf-8") as handle:
                for line in handle:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    if line.startswith("Annotation\t"):
                        continue
                    out.write(f"{source}\t{line}\n")


if __name__ == "__main__":
    main()
