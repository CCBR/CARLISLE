#!/usr/bin/env python3
"""Aggregate treatment-level ROSE outputs across treatments."""

from __future__ import annotations

import argparse
import os


def parse_input_spec(spec: str):
    if "::" not in spec:
        raise ValueError(f"Invalid --input '{spec}', expected label::path")
    label, path = spec.split("::", 1)
    return label, path


def parse_label(label: str):
    parts = label.split("|")
    if len(parts) != 5:
        raise ValueError(
            f"Invalid aggregate label '{label}', expected peak_caller|control_mode|treatment_sample|dup_status|stitch_distance"
        )
    return parts


def aggregate(input_specs, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as out:
        out.write("peak_caller\ttreatment_sample\tcontrol_mode\tdup_status\tstitch_distance\trow\n")
        for spec in input_specs:
            source, path = parse_input_spec(spec)
            peak_caller, control_mode, treatment_sample, dup_status, stitch_distance = parse_label(source)
            if not os.path.exists(path):
                continue
            with open(path, "r", encoding="utf-8") as handle:
                for line in handle:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith("Less than 5 usable peaks detected"):
                        continue
                    out.write(
                        f"{peak_caller}\t{treatment_sample}\t{control_mode}\t{dup_status}\t{stitch_distance}\t{line}\n"
                    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate treatment-level ROSE mapping outputs")
    parser.add_argument("--enhancer-to-gene", action="append", default=[])
    parser.add_argument("--gene-to-enhancer", action="append", default=[])
    parser.add_argument("--output-enhancer-to-gene", required=True)
    parser.add_argument("--output-gene-to-enhancer", required=True)
    args = parser.parse_args()

    aggregate(args.enhancer_to_gene, args.output_enhancer_to_gene)
    aggregate(args.gene_to_enhancer, args.output_gene_to_enhancer)


if __name__ == "__main__":
    main()
