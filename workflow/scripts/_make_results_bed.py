#!/usr/bin/env python3
"""
Create a colored BED file from differential analysis results.

This script processes differential analysis results and generates a BED file with
color-coded peaks based on log2 fold change and adjusted p-value thresholds. The
colors indicate the significance and direction of differential binding/expression.

The script applies the following color scheme:
- Gray (160,160,160): Non-significant peaks
- Light green (229,255,204): Significantly downregulated (between elbow and log2FC cutoff)
- Bright green (128,255,0): Highly downregulated (below log2FC cutoff)
- Orange (255,153,51): Significantly upregulated (between elbow and log2FC cutoff)
- Red (255,51,51): Highly upregulated (above log2FC cutoff)

The strand field is set based on the direction of fold change:
- '+' for upregulated (log2FoldChange >= 0)
- '-' for downregulated (log2FoldChange < 0)

Command-line Arguments:
    --results (str): Path to the input results file (TSV format) containing
                    differential analysis results with columns: seqnames, start,
                    end, log2FoldChange, and padj.
    --fdr_cutoff (float): False Discovery Rate (FDR) threshold for significance.
                         Peaks with padj < fdr_cutoff are considered significant.
    --log2FC_cutoff (float): Log2 fold change threshold for high significance.
                            Peaks exceeding this threshold are considered highly
                            differentially expressed.
    --elbowyaml (str): Path to YAML file containing elbow point limits with
                      'low_limit' (downregulation elbow) and 'up_limit'
                      (upregulation elbow) keys.
    --bed (str): Path to the output BED file with color annotations.

Output:
    A sorted BED file with 9 columns:
    1. seqnames (chromosome)
    2. start (start position)
    3. end (end position)
    4. name (always '.')
    5. score (always 0)
    6. strand ('+' or '-' based on fold change direction)
    7. seven (always 0)
    8. eight (always 0)
    9. color (RGB color based on significance and fold change)

Note:
    - The script uses a temporary file with a UUID-based name to avoid conflicts
    - The final output is sorted by chromosome and start position
    - Sorting uses 40GB of memory and /dev/shm for temporary storage
    - The temporary file is automatically cleaned up after processing
"""

import argparse, pandas, os, yaml, subprocess, uuid
from operator import index
from email import header

randbed = str(uuid.uuid4()) + ".bed"

parser = argparse.ArgumentParser(description="create bed")
parser.add_argument("--results", required=True, type=str, help="results bed file")
parser.add_argument("--fdr_cutoff", required=True, type=float, help="FDR cutoff")
parser.add_argument("--log2FC_cutoff", required=True, type=float, help="log2FC cutoff")
parser.add_argument("--elbowyaml", required=True, type=str, help="elbow limits")
parser.add_argument("--bed", required=True, type=str, help="output bed with colors")
args = parser.parse_args()

df = pandas.read_csv(
    args.results,
    sep="\t",
    header=0,
    usecols=["seqnames", "start", "end", "log2FoldChange", "padj"],
)
df["color"] = ["160,160,160"] * len(df.index)
df["name"] = ["."] * len(df.index)
df["strand"] = ["+"] * len(df.index)
df["score"] = [0] * len(df.index)
df["seven"] = [0] * len(df.index)
df["eight"] = [0] * len(df.index)

df.loc[df["log2FoldChange"] < 0, "strand"] = "-"
with open(args.elbowyaml) as f:
    elbowdf = yaml.safe_load(f)
elbowdown = elbowdf["low_limit"]
elbowup = elbowdf["up_limit"]
if abs(elbowdown) < args.log2FC_cutoff:
    df.loc[
        (df["padj"] < args.fdr_cutoff) & (df["log2FoldChange"] < elbowdown), "color"
    ] = "229,255,204"
df.loc[
    (df["padj"] < args.fdr_cutoff) & (df["log2FoldChange"] < args.log2FC_cutoff),
    "color",
] = "128,255,0"
if elbowup < args.log2FC_cutoff:
    df.loc[
        (df["padj"] < args.fdr_cutoff) & (df["log2FoldChange"] > elbowup), "color"
    ] = "255,153,51"
df.loc[
    (df["padj"] < args.fdr_cutoff) & (df["log2FoldChange"] > args.log2FC_cutoff),
    "color",
] = "255,51,51"
df = df.reindex(
    columns=[
        "seqnames",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "seven",
        "eight",
        "color",
    ]
)

df.to_csv(randbed, sep="\t", header=False, index=False)

cmd = "sort -S 40G -T /dev/shm -k1,1 -k2,2n " + randbed + " > " + args.bed
s = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)

os.remove(randbed)
