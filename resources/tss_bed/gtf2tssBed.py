#!/usr/bin/env python3
import argparse
import os

def get_gene_name(s):
        s=list(map(lambda x:x.replace(";","").replace("\"",""), s.strip().split(" ")))
        for i,v in enumerate(s):
                if v == "gene_name":
                        gene_name = s[i+1]
                        return gene_name
                elif v == "gene":
                        gene_name = s[i+1]
                        return gene_name
        else:
                return "Unknown"

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", help="GTF input file", required=True)
parser.add_argument("--bed", help="TSS BED output file", required=True)
args = parser.parse_args()

def get_gene_start(start,end,strand):
        if strand == "+":
                return int(start)-1
        else:
                return int(end)-1

if not os.path.exists(args.gtf):
        print("%s : GTF File is missing!"%(args.gtf))
        exit()

gtf = open(args.gtf,'r')
bed = open(args.bed,'w')
for l in gtf.readlines():
        if l.startswith("#"): continue
        l = l.split("\t")
        if l[2] == "gene":
                gene_name = get_gene_name(l[8])
                chrom = l[0]
                start = l[3]
                end = l[4]
                strand = l[6]
                gene_start = get_gene_start(start,end,strand)
                bed_line = [chrom, str(gene_start), str(gene_start+1),gene_name,str(0),strand]
                bed.write("%s\n"%("\t".join(bed_line)))
gtf.close()
bed.close()
