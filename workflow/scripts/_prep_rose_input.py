#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Prepare and validate ROSE inputs from peak calls (Python 2.7 compatible).

This script covers ROSE preflight items:
1) Peak input normalization (MACS2/SEACR/GoPeaks/BED -> BED6)
2) Alignment input checks (treatment/control BAM + index)
3) Reference asset checks (TSS/optional blacklist)
4) ROSE parameter validation (stitch/tss distance, labels)
5) Software/runtime checks (bedtools/samtools/ROSE_main.py)
6) Pre-ROSE filtering (blacklist removal + TSS exclusion + stitching)
7) Resource estimate summary
8) Expected-output manifest (and optional ROSE execution)
"""

from __future__ import print_function

import argparse
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile


def run_cmd(cmd, fail_msg, capture=True, cwd=None):
    try:
        if capture:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            out, err = proc.communicate()
            rc = proc.returncode
        else:
            proc = subprocess.Popen(cmd, cwd=cwd)
            rc = proc.wait()
            out, err = "", ""
        if rc != 0:
            raise RuntimeError(
                "%s\nCommand: %s\n%s" % (
                    fail_msg,
                    " ".join(cmd),
                    (err or "").strip(),
                )
            )
        return out, err
    except OSError as exc:
        raise RuntimeError("%s\nCommand: %s\n%s" % (fail_msg, " ".join(cmd), str(exc)))


def ensure_file(path, label):
    if not os.path.exists(path):
        raise IOError("%s does not exist: %s" % (label, path))
    if os.path.isdir(path):
        raise IOError("%s points to a directory, expected file: %s" % (label, path))


def ensure_executable(name_or_path, label):
    if os.path.sep in name_or_path:
        if os.path.exists(name_or_path):
            return name_or_path
        raise IOError("%s not found: %s" % (label, name_or_path))
    resolved = shutil.which(name_or_path) if hasattr(shutil, "which") else None
    if not resolved:
        # py2 fallback
        for p in os.environ.get("PATH", "").split(os.pathsep):
            candidate = os.path.join(p, name_or_path)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                resolved = candidate
                break
    if not resolved:
        raise IOError("%s not found in PATH: %s" % (label, name_or_path))
    return resolved


def parse_number(value, default=0.0):
    try:
        val = float(value)
        if math.isfinite(val):
            return val
    except Exception:
        pass
    return default


def sanitize_sample_id(sample_id):
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", sample_id.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    if not safe:
        raise ValueError("sample_id is empty after sanitization")
    return safe


def detect_peak_format(path):
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 10:
                return "macs2_narrowpeak"
            if len(fields) >= 6:
                if ":" in fields[5] and "-" in fields[5]:
                    return "seacr"
                return "bed"
            if len(fields) >= 3:
                return "gopeaks"
    return "bed"


def load_peaks(path, peak_format, sample_id):
    rows = []
    sample = sanitize_sample_id(sample_id)
    idx = 0
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                start = int(float(fields[1]))
                end = int(float(fields[2]))
            except ValueError:
                continue
            if start < 0 or end <= start:
                continue

            idx += 1
            name = "%s_%06d" % (sample, idx)
            score = 0.0
            strand = "."

            if peak_format == "seacr" and len(fields) >= 4:
                score = parse_number(fields[3], default=0.0)
            if len(fields) >= 5:
                score = parse_number(fields[4], default=score)
            if len(fields) >= 6 and fields[5] in ("+", "-"):
                strand = fields[5]

            rows.append((chrom, start, end, name, score, strand))

    return rows


def write_bed6(rows, out_bed):
    rows_sorted = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    out_dir = os.path.dirname(out_bed)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(out_bed, "w") as handle:
        for r in rows_sorted:
            handle.write("%s\t%d\t%d\t%s\t%.6f\t%s\n" % (r[0], r[1], r[2], r[3], r[4], r[5]))


def maybe_decompress_gz(in_path, workdir):
    if not in_path.endswith(".gz"):
        return in_path
    out_path = os.path.join(workdir, os.path.basename(in_path[:-3]))
    with gzip.open(in_path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def filter_with_bedtools(bedtools_bin, in_bed, exclude_bed, out_bed, reason):
    cmd = [bedtools_bin, "intersect", "-a", in_bed, "-b", exclude_bed, "-v"]
    out, _ = run_cmd(cmd, fail_msg="bedtools intersect failed during %s" % reason, capture=True)
    with open(out_bed, "w") as handle:
        handle.write(out.decode("utf-8") if isinstance(out, bytes) else out)


def stitch_with_bedtools(bedtools_bin, in_bed, stitch_distance, out_bed):
    sort_cmd = [bedtools_bin, "sort", "-i", in_bed]
    out, _ = run_cmd(sort_cmd, fail_msg="bedtools sort failed for stitched ROSE input", capture=True)
    if isinstance(out, bytes):
        out = out.decode("utf-8")
    fd, tmp_path = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        with open(tmp_path, "w") as tmp:
            tmp.write(out)
        merge_cmd = [
            bedtools_bin,
            "merge",
            "-i",
            tmp_path,
            "-d",
            str(stitch_distance),
            "-c",
            "4,5,6",
            "-o",
            "distinct,sum,distinct",
        ]
        merged, _ = run_cmd(merge_cmd, fail_msg="bedtools merge failed while stitching peaks", capture=True)
        if isinstance(merged, bytes):
            merged = merged.decode("utf-8")
        with open(out_bed, "w") as handle:
            handle.write(merged)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def convert_stitched_bed_to_rose_gff(in_bed, out_gff):
    """
    Convert stitched BED6 to ROSE-friendly GFF-like table using:
    chrom, name, ., start(1-based), end, ., ., ., ., ., ., name
    """
    with open(in_bed, "r") as src, open(out_gff, "w") as dst:
        rownum = 0
        for raw in src:
            line = raw.strip()
            if not line:
                continue
            rownum += 1
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start0 = int(float(fields[1]))
            end = int(float(fields[2]))
            name = fields[3] if len(fields) >= 4 and fields[3] else "region%d" % rownum
            start1 = start0 + 1
            dst.write(
                "%s\t%s\t.\t%d\t%d\t.\t.\t.\t.\t.\t.\t%s\n"
                % (chrom, name, start1, end, name)
            )


def count_nonempty_lines(path):
    n = 0
    with open(path, "r") as handle:
        for line in handle:
            if line.strip():
                n += 1
    return n


def estimate_resources(stitched_peak_count):
    if stitched_peak_count < 10000:
        return {"threads": 1, "mem": "16g", "time": "02:00:00"}
    if stitched_peak_count < 50000:
        return {"threads": 2, "mem": "32g", "time": "04:00:00"}
    if stitched_peak_count < 200000:
        return {"threads": 4, "mem": "64g", "time": "08:00:00"}
    return {"threads": 8, "mem": "120g", "time": "16:00:00"}


def expected_rose_outputs(output_dir, sample_id):
    prefix = sanitize_sample_id(sample_id)
    return [
        os.path.join(output_dir, "%s_AllStitched.table.txt" % prefix),
        os.path.join(output_dir, "%s_AllEnhancers.table.txt" % prefix),
        os.path.join(output_dir, "%s_SuperEnhancers.table.txt" % prefix),
        os.path.join(output_dir, "%s_Enhancer_to_Gene.table.txt" % prefix),
    ]


def write_reports(output_dir, report):
    report_json = os.path.join(output_dir, "rose_preflight_report.json")
    report_txt = os.path.join(output_dir, "rose_preflight_report.txt")
    with open(report_json, "w") as outj:
        json.dump(report, outj, indent=2, sort_keys=True)

    with open(report_txt, "w") as outt:
        outt.write("ROSE preflight summary\n")
        outt.write("======================\n")
        outt.write("Sample: %s\n" % report["sample_id"])
        outt.write("Genome: %s\n" % report["genome"])
        outt.write("Peak format: %s\n" % report["peak_format"])
        outt.write("Input peaks: %d\n" % report["counts"]["input_bed6"])
        outt.write("After TSS exclusion: %d\n" % report["counts"]["after_tss_exclusion"])
        outt.write("Stitched peaks: %d\n" % report["counts"]["stitched"])
        outt.write("Prepared stitched BED: %s\n" % report["prepared_stitched_bed"])
        outt.write("Prepared stitched GFF: %s\n" % report["prepared_stitched_gff"])
        outt.write("ROSE runnable (>= min_peaks=%d): %s\n" % (report["min_peaks"], str(report["run_rose_allowed"])))
        outt.write(
            "Resource estimate: threads=%s, mem=%s, time=%s\n"
            % (
                report["resource_estimate"]["threads"],
                report["resource_estimate"]["mem"],
                report["resource_estimate"]["time"],
            )
        )
        outt.write("Expected outputs:\n")
        for item in report["expected_outputs"]:
            outt.write("  - %s\n" % item)
        if report.get("notes"):
            outt.write("Notes:\n")
            for n in report["notes"]:
                outt.write("  - %s\n" % n)
    return report_json, report_txt


def validate_bam_and_index(bam_path, label):
    ensure_file(bam_path, label)
    bai1 = bam_path + ".bai"
    bai2 = re.sub(r"\.bam$", ".bai", bam_path)
    if not (os.path.exists(bai1) or os.path.exists(bai2)):
        raise IOError("%s BAM index missing (.bai): %s" % (label, bam_path))


def main():
    parser = argparse.ArgumentParser(
        description="Prepare ROSE-compatible stitched BED from MACS2/SEACR/GoPeaks peaks and run preflight checks."
    )
    parser.add_argument("--peak-file", required=True, help="Input peak file (MACS2 narrow/broadPeak, SEACR, GoPeaks, or BED).")
    parser.add_argument(
        "--peak-format",
        default="auto",
        choices=["auto", "macs2_narrowpeak", "macs2_broadpeak", "seacr", "gopeaks", "bed"],
        help="Explicit peak format; default auto-detect.",
    )
    parser.add_argument("--sample-id", required=True, help="Sample label used for ROSE region IDs and expected outputs.")
    parser.add_argument("--treatment-bam", required=True, help="Treatment BAM used by ROSE (-r).")
    parser.add_argument("--control-bam", default=None, help="Optional control BAM used by ROSE (-r).")
    parser.add_argument("--genome", required=True, help="Genome label (e.g., hg38, mm10) for reporting.")
    parser.add_argument("--tss-bed", required=True, help="TSS BED file used for promoter exclusion (can be .gz).")
    parser.add_argument("--blacklist-bed", default=None, help="Optional blacklist BED to remove problematic loci.")
    parser.add_argument("--stitch-distance", type=int, required=True, help="ROSE stitching distance in bp.")
    parser.add_argument("--tss-distance", type=int, required=True, help="ROSE TSS exclusion distance in bp.")
    parser.add_argument("--output-dir", required=True, help="Output directory for prepared files and reports.")
    parser.add_argument("--prepared-bed-name", default="rose_input.prepared.stitched.bed", help="Prepared stitched BED filename.")
    parser.add_argument("--prepared-gff-name", default="rose_input.prepared.stitched.gff", help="Prepared stitched GFF filename.")
    parser.add_argument("--min-peaks", type=int, default=5, help="Minimum stitched peak count required to run ROSE.")
    parser.add_argument("--rose-main", default="ROSE_main.py", help="ROSE executable name/path.")
    parser.add_argument("--rose-root", default="/opt/ROSE", help="ROSE installation root; ROSE_main.py is run with cwd here.")
    parser.add_argument("--rose-python", default="/opt/conda/envs/rose/bin/python", help="Python executable used to run ROSE_main.py.")
    parser.add_argument("--run-rose", action="store_true", help="Run ROSE_main.py after preflight and BED preparation.")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate BED files in output dir.")
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    intermediate_dir = os.path.join(output_dir, "rose_prep_intermediate")
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    peak_file = os.path.abspath(args.peak_file)
    tss_bed = os.path.abspath(args.tss_bed)
    treatment_bam = os.path.abspath(args.treatment_bam)
    control_bam = os.path.abspath(args.control_bam) if args.control_bam else None
    blacklist_bed = None
    if args.blacklist_bed and str(args.blacklist_bed).strip().lower() not in ("", "none", "na", "null"):
        blacklist_bed = os.path.abspath(args.blacklist_bed)

    ensure_file(peak_file, "peak file")
    ensure_file(tss_bed, "tss bed")
    if blacklist_bed:
        ensure_file(blacklist_bed, "blacklist bed")

    validate_bam_and_index(treatment_bam, "treatment")
    if control_bam:
        validate_bam_and_index(control_bam, "control")

    if args.stitch_distance <= 0 or args.tss_distance < 0:
        raise ValueError("stitch_distance must be > 0 and tss_distance must be >= 0")
    if args.min_peaks < 1:
        raise ValueError("min_peaks must be >= 1")

    bedtools_bin = ensure_executable("bedtools", "bedtools")
    samtools_bin = ensure_executable("samtools", "samtools")
    rose_py = ensure_executable(args.rose_python, "ROSE python") if args.run_rose else args.rose_python
    rose_root = os.path.abspath(args.rose_root)
    rose_main_path = args.rose_main
    if args.run_rose:
        if os.path.sep not in args.rose_main:
            rose_main_path = os.path.join(rose_root, args.rose_main)
        else:
            rose_main_path = os.path.abspath(args.rose_main)
        ensure_file(rose_main_path, "ROSE_main.py")
        if not os.path.isdir(rose_root):
            raise IOError("rose-root does not exist or is not a directory: %s" % rose_root)

    bedtools_ver, _ = run_cmd([bedtools_bin, "--version"], fail_msg="Unable to query bedtools version", capture=True)
    samtools_ver, _ = run_cmd([samtools_bin, "--version"], fail_msg="Unable to query samtools version", capture=True)
    if isinstance(bedtools_ver, bytes):
        bedtools_ver = bedtools_ver.decode("utf-8")
    if isinstance(samtools_ver, bytes):
        samtools_ver = samtools_ver.decode("utf-8")
    bedtools_ver = bedtools_ver.strip()
    samtools_ver = samtools_ver.splitlines()[0].strip() if samtools_ver else ""

    peak_format = args.peak_format
    if peak_format == "auto":
        peak_format = detect_peak_format(peak_file)

    raw_rows = load_peaks(peak_file, peak_format=peak_format, sample_id=args.sample_id)
    bed6_raw = os.path.join(intermediate_dir, "01_input_as_bed6.bed")
    write_bed6(raw_rows, bed6_raw)

    raw_n = len(raw_rows)
    stitched_bed = os.path.join(output_dir, args.prepared_bed_name)
    stitched_gff = os.path.join(output_dir, args.prepared_gff_name)

    # Early exit for empty/tiny peak sets before expensive filtering/stitching.
    if raw_n < args.min_peaks:
        open(stitched_bed, "w").close()
        open(stitched_gff, "w").close()
        report = {
            "genome": args.genome,
            "sample_id": sanitize_sample_id(args.sample_id),
            "peak_file": peak_file,
            "peak_format": peak_format,
            "treatment_bam": treatment_bam,
            "control_bam": control_bam,
            "tss_bed": tss_bed,
            "blacklist_bed": blacklist_bed,
            "stitch_distance": args.stitch_distance,
            "tss_distance": args.tss_distance,
            "min_peaks": args.min_peaks,
            "counts": {
                "input_bed6": raw_n,
                "after_blacklist": raw_n,
                "after_tss_exclusion": 0,
                "stitched": 0,
            },
            "software": {
                "bedtools": bedtools_ver,
                "samtools": samtools_ver,
                "rose_main": args.rose_main,
            },
            "resource_estimate": estimate_resources(0),
            "expected_outputs": expected_rose_outputs(output_dir, args.sample_id),
            "prepared_stitched_bed": stitched_bed,
            "prepared_stitched_gff": stitched_gff,
            "run_rose_requested": args.run_rose,
            "run_rose_allowed": False,
            "rose_command": None,
            "notes": [
                "Input peaks below min_peaks threshold; skipping filtering/stitching and ROSE execution.",
            ],
        }
        report_json, report_txt = write_reports(output_dir, report)
        if not args.keep_intermediate and os.path.isdir(intermediate_dir):
            for p in os.listdir(intermediate_dir):
                fp = os.path.join(intermediate_dir, p)
                if os.path.isfile(fp):
                    os.remove(fp)
            try:
                os.rmdir(intermediate_dir)
            except OSError:
                pass
        print("[OK] Peak count below threshold (%d < %d); wrote empty prepared files." % (raw_n, args.min_peaks))
        print("[OK] Prepared ROSE stitched BED: %s" % stitched_bed)
        print("[OK] Prepared ROSE stitched GFF: %s" % stitched_gff)
        print("[OK] Wrote report JSON: %s" % report_json)
        print("[OK] Wrote report TXT:  %s" % report_txt)
        return 0

    stage_bed = bed6_raw
    if blacklist_bed:
        no_blacklist = os.path.join(intermediate_dir, "02_no_blacklist.bed")
        filter_with_bedtools(bedtools_bin, stage_bed, blacklist_bed, no_blacklist, "blacklist filtering")
        stage_bed = no_blacklist

    tss_plain = maybe_decompress_gz(tss_bed, intermediate_dir)
    no_tss = os.path.join(intermediate_dir, "03_no_tss_overlap.bed")
    filter_with_bedtools(bedtools_bin, stage_bed, tss_plain, no_tss, "TSS exclusion filtering")

    stitch_with_bedtools(bedtools_bin, no_tss, args.stitch_distance, stitched_bed)
    convert_stitched_bed_to_rose_gff(stitched_bed, stitched_gff)

    raw_n = count_nonempty_lines(bed6_raw)
    no_blacklist_n = count_nonempty_lines(stage_bed) if os.path.exists(stage_bed) else raw_n
    no_tss_n = count_nonempty_lines(no_tss)
    stitched_n = count_nonempty_lines(stitched_bed)

    run_rose_allowed = stitched_n >= args.min_peaks
    if args.run_rose and not run_rose_allowed:
        raise RuntimeError(
            "Stitched peaks below threshold (%d < %d); refusing to run ROSE. Prepared BED: %s"
            % (stitched_n, args.min_peaks, stitched_bed)
        )

    rose_cmd = None
    if args.run_rose:
        rose_cmd = [
            rose_py,
            rose_main_path,
            "-g",
            args.genome,
            "-i",
            stitched_gff,
            "-r",
            treatment_bam,
            "-s",
            str(args.stitch_distance),
            "-t",
            str(args.tss_distance),
            "-o",
            output_dir,
        ]
        if control_bam:
            rose_cmd.extend(["-c", control_bam])
        run_cmd(rose_cmd, fail_msg="ROSE_main.py failed", capture=True, cwd=rose_root)

    report = {
        "genome": args.genome,
        "sample_id": sanitize_sample_id(args.sample_id),
        "peak_file": peak_file,
        "peak_format": peak_format,
        "treatment_bam": treatment_bam,
        "control_bam": control_bam,
        "tss_bed": tss_bed,
        "blacklist_bed": blacklist_bed,
        "stitch_distance": args.stitch_distance,
        "tss_distance": args.tss_distance,
        "min_peaks": args.min_peaks,
        "counts": {
            "input_bed6": raw_n,
            "after_blacklist": no_blacklist_n,
            "after_tss_exclusion": no_tss_n,
            "stitched": stitched_n,
        },
        "software": {
            "bedtools": bedtools_ver,
            "samtools": samtools_ver,
            "rose_main": args.rose_main,
        },
        "resource_estimate": estimate_resources(stitched_n),
        "expected_outputs": expected_rose_outputs(output_dir, args.sample_id),
        "prepared_stitched_bed": stitched_bed,
        "prepared_stitched_gff": stitched_gff,
        "run_rose_requested": args.run_rose,
        "run_rose_allowed": run_rose_allowed,
        "rose_command": rose_cmd,
    }
    report_json, report_txt = write_reports(output_dir, report)

    if not args.keep_intermediate and os.path.isdir(intermediate_dir):
        for p in os.listdir(intermediate_dir):
            fp = os.path.join(intermediate_dir, p)
            if os.path.isfile(fp):
                os.remove(fp)
        try:
            os.rmdir(intermediate_dir)
        except OSError:
            pass

    print("[OK] Prepared ROSE stitched BED: %s" % stitched_bed)
    print("[OK] Prepared ROSE stitched GFF: %s" % stitched_gff)
    print("[OK] Wrote report JSON: %s" % report_json)
    print("[OK] Wrote report TXT:  %s" % report_txt)
    if args.run_rose:
        print("[OK] ROSE execution finished.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("[ERROR] %s" % str(exc), file=sys.stderr)
        sys.exit(1)
