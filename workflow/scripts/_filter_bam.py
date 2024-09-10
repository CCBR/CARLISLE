#!/usr/bin/env python3
# Author: Vishal Koparde, PhD
# Date: Jan 2022
# This script take raw alignments, presorted and indexed and filtered them for reads
# aligning as proper pairs with a set fragment size.
# inputs:
# @inputBAM: raw BAM file, sorted and indexed
# @flagmentlength: default 1000, readpairs with fragment length larger than this integer are discarded
# outputs:
# @outputBAM: filtered output BAM file
import pysam
import argparse

parser = argparse.ArgumentParser(description="Filter BAM by readids")
parser.add_argument("--inputBAM", type=str, required=True, help="input BAM file")
parser.add_argument(
    "--outputBAM", type=str, required=True, help="filtered output BAM file"
)
parser.add_argument(
    "--fragmentlength",
    type=int,
    required=False,
    default=1000,
    help="discard flagment lengths larger than this integer",
)
parser.add_argument(
    "--removemarkedduplicates",
    action="store_true",
    help="removed marked and optical duplicates",
)
args = parser.parse_args()

inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)

count = 0
for read in inBAM.fetch():
    count += 1
    if count % 1000000 == 0:
        print("%d reads read!" % (count))
    if not read.is_proper_pair or read.is_unmapped:
        continue
    if read.template_length > args.fragmentlength:
        continue
    if args.removemarkedduplicates and read.is_duplicate:
        continue
    outBAM.write(read)
inBAM.close()
outBAM.close()
