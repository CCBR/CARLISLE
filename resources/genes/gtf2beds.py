#!/usr/bin/env python3
import argparse
import os


def parse_attributes(attr_str):
    """
    Parse GTF attribute column into a dict.
    Example: 'gene_id "X"; transcript_id "Y"; gene_name "Z";'
    """
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.replace('"', "").split()
        if len(parts) >= 2:
            key = parts[0]
            value = parts[1]
            attrs[key] = value
    return attrs


def get_tss(min_start, max_end, strand):
    """
    Given aggregated gene min_start/max_end (1-based GTF coords),
    return TSS (1-based).
    """
    return min_start if strand == "+" else max_end


def load_genome_sizes(genome_sizes_file):
    """
    Load genome sizes into dict: chrom -> length (int).
    """
    gs = {}
    with open(genome_sizes_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            length = int(parts[1])
            gs[chrom] = length
    return gs


def load_bed_intervals(bed_file, genome_sizes):
    """
    Load BED intervals into dict: chrom -> list of (start, end).
    Only keeps chromosomes present in genome_sizes.
    """
    intervals = {chrom: [] for chrom in genome_sizes.keys()}
    if bed_file is None:
        return intervals

    with open(bed_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in genome_sizes:
                continue
            start = int(parts[1])
            end = int(parts[2])
            if end <= start:
                continue
            intervals[chrom].append((start, end))
    # Merge per chrom for efficiency
    for chrom in intervals:
        intervals[chrom] = merge_intervals(intervals[chrom])
    return intervals


def merge_intervals(intervals):
    """
    Merge overlapping or adjacent intervals.
    intervals: list of (start, end) with 0-based half-open coordinates.
    Returns a new list of merged intervals.
    """
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]

    for curr_start, curr_end in intervals[1:]:
        last_start, last_end = merged[-1]
        if curr_start <= last_end:  # overlap or contiguous
            merged[-1] = (last_start, max(last_end, curr_end))
        else:
            merged.append((curr_start, curr_end))

    return merged


def subtract_intervals(base_intervals, subtract_intervals):
    """
    Subtract a set of intervals from base_intervals.
    All intervals are 0-based half-open and assumed sorted/merged per set.
    base_intervals: list[(start, end)]
    subtract_intervals: list[(start, end)]
    Returns new list of base intervals with subtract_intervals removed.
    """
    if not base_intervals or not subtract_intervals:
        return base_intervals[:]

    result = []
    j = 0
    n_sub = len(subtract_intervals)

    for b_start, b_end in base_intervals:
        curr_start = b_start
        while j < n_sub and subtract_intervals[j][1] <= curr_start:
            # sub interval ends before base interval starts
            j += 1
        k = j
        while k < n_sub and subtract_intervals[k][0] < b_end:
            s_start, s_end = subtract_intervals[k]
            # no overlap ahead
            if s_start > curr_start:
                result.append((curr_start, min(s_start, b_end)))
            curr_start = max(curr_start, s_end)
            if curr_start >= b_end:
                break
            k += 1
        if curr_start < b_end:
            result.append((curr_start, b_end))

    # Filter any empties
    result = [(s, e) for s, e in result if e > s]
    return result

parser = argparse.ArgumentParser(
    description="Generate TSS, promoter, genebody, and intergenic BED files "
                "from a GTF file, optionally excluding blacklist and ambiguous regions."
)

# =========================
# INPUTS
# =========================
input_group = parser.add_argument_group("INPUT FILES")

input_group.add_argument(
    "--gtf",
    required=True,
    help="INPUT: GTF annotation file"
)

input_group.add_argument(
    "--genome_sizes",
    required=True,
    help="INPUT: Genome sizes file (chrom\\tlength)"
)

input_group.add_argument(
    "--blacklist_bed",
    required=False,
    help="INPUT: BED file of blacklist regions to exclude from intergenic regions"
)

input_group.add_argument(
    "--ambiguous_bed",
    required=False,
    help="INPUT: BED file of ambiguous-base regions (e.g., Ns) to exclude from intergenic"
)

input_group.add_argument(
    "--promoter_window",
    type=int,
    default=1000,
    help="INPUT: Half-window size around TSS (default = 1000 for +/-1kb)"
)

# =========================
# OUTPUTS
# =========================
output_group = parser.add_argument_group("OUTPUT FILES")

output_group.add_argument(
    "--tss_bed",
    required=True,
    help="OUTPUT: BED file for TSS (1bp)"
)

output_group.add_argument(
    "--promoter_bed",
    required=True,
    help="OUTPUT: BED file for promoter regions (+/- promoter_window)"
)

output_group.add_argument(
    "--genebody_bed",
    required=True,
    help="OUTPUT: BED file for longest gene bodies (merged isoforms)"
)

output_group.add_argument(
    "--intergenic_bed",
    required=True,
    help="OUTPUT: BED file for intergenic regions "
         "(genome - genebody - blacklist - ambiguous)"
)

args = parser.parse_args()
# ---- Sanity check inputs ----
if not os.path.exists(args.gtf):
    print(f"{args.gtf} : GTF file is missing!")
    exit(1)

if not os.path.exists(args.genome_sizes):
    print(f"{args.genome_sizes} : Genome sizes file is missing!")
    exit(1)

if args.blacklist_bed and not os.path.exists(args.blacklist_bed):
    print(f"{args.blacklist_bed} : Blacklist BED file is missing!")
    exit(1)

if args.ambiguous_bed and not os.path.exists(args.ambiguous_bed):
    print(f"{args.ambiguous_bed} : Ambiguous BED file is missing!")
    exit(1)

# ---- Load genome sizes ----
genome_sizes = load_genome_sizes(args.genome_sizes)

# ---- Load blacklist and ambiguous intervals ----
blacklist_intervals = load_bed_intervals(args.blacklist_bed, genome_sizes)
ambiguous_intervals = load_bed_intervals(args.ambiguous_bed, genome_sizes)

# ---- First pass: aggregate gene spans (longest isoform) ----
genes = {}  # gene_id -> dict

with open(args.gtf, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        chrom = cols[0]
        feature = cols[2]
        start = int(cols[3])
        end = int(cols[4])
        strand = cols[6]
        attr_str = cols[8]

        attrs = parse_attributes(attr_str)
        gene_id = attrs.get("gene_id")
        if gene_id is None:
            continue  # skip entries without gene_id

        gene_name = attrs.get("gene_name", gene_id)

        # Features contributing to overall gene span:
        if feature not in {"gene", "transcript", "mRNA", "exon", "CDS"}:
            continue

        if gene_id not in genes:
            genes[gene_id] = {
                "chrom": chrom,
                "strand": strand,
                "min_start": start,
                "max_end": end,
                "gene_name": gene_name,
            }
        else:
            g = genes[gene_id]
            if start < g["min_start"]:
                g["min_start"] = start
            if end > g["max_end"]:
                g["max_end"] = end

# ---- Second pass: write TSS, promoter, gene body, and collect intervals per chrom ----
w = args.promoter_window

# For building intergenic: chrom -> list of (start, end) for genebody
genebody_intervals = {chrom: [] for chrom in genome_sizes.keys()}

with open(args.promoter_bed, "w") as promoter_out, \
     open(args.tss_bed, "w") as tss_out, \
     open(args.genebody_bed, "w") as genebody_out:

    for gene_id, info in genes.items():
        chrom = info["chrom"]
        strand = info["strand"]
        min_start = info["min_start"]      # 1-based
        max_end = info["max_end"]          # 1-based inclusive
        gene_name = info["gene_name"]

        # Skip genes on chromosomes not present in genome_sizes
        if chrom not in genome_sizes:
            continue

        # --- TSS (1-based -> 0-based BED) ---
        tss_1based = get_tss(min_start, max_end, strand)
        tss_0based = tss_1based - 1

        tss_out.write(
            "\t".join([
                chrom,
                str(tss_0based),
                str(tss_0based + 1),
                gene_name,
                "0",
                strand
            ]) + "\n"
        )

        # --- Promoter: +/- w around TSS (0-based, half-open) ---
        promoter_start = max(0, tss_0based - w)
        promoter_end = tss_0based + w  # BED end is exclusive

        promoter_out.write(
            "\t".join([
                chrom,
                str(promoter_start),
                str(promoter_end),
                gene_name,
                "0",
                strand
            ]) + "\n"
        )

        # --- Gene body: full span of gene (longest isoform) ---
        # Convert 1-based [min_start, max_end] to 0-based BED [start, end)
        genebody_start = min_start - 1
        genebody_end = max_end  # BED end is exclusive

        if genebody_start < 0:
            genebody_start = 0

        genebody_out.write(
            "\t".join([
                chrom,
                str(genebody_start),
                str(genebody_end),
                gene_name,
                "0",
                strand
            ]) + "\n"
        )

        genebody_intervals[chrom].append((genebody_start, genebody_end))

# Merge gene body intervals per chromosome
for chrom in genebody_intervals:
    genebody_intervals[chrom] = merge_intervals(genebody_intervals[chrom])

# ---- Third pass: compute intergenic regions as complement of genebody, then subtract blacklist & ambiguous ----
with open(args.intergenic_bed, "w") as intergenic_out:
    for chrom, chrom_len in genome_sizes.items():
        # Start with complement of genebody
        gb = genebody_intervals.get(chrom, [])
        intergenic = []

        if not gb:
            # Whole chromosome is initially intergenic
            intergenic = [(0, chrom_len)]
        else:
            prev_end = 0
            for start, end in gb:
                if start > prev_end:
                    intergenic.append((prev_end, start))
                prev_end = max(prev_end, end)
            if prev_end < chrom_len:
                intergenic.append((prev_end, chrom_len))

        # Subtract blacklist and ambiguous intervals
        bl = blacklist_intervals.get(chrom, [])
        amb = ambiguous_intervals.get(chrom, [])

        if bl:
            intergenic = subtract_intervals(intergenic, bl)
        if amb:
            intergenic = subtract_intervals(intergenic, amb)

        # Write out final intergenic intervals
        for start, end in intergenic:
            if end <= start:
                continue
            intergenic_out.write(
                "\t".join([
                    chrom,
                    str(start),
                    str(end),
                    "intergenic",
                    "0",
                    "."
                ]) + "\n"
            )

