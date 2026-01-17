#!/usr/bin/env bash
set -euo pipefail

OUT_ROOT="/data/CCBR_Pipeliner/Pipelines/CARLISLE/vnk_dev/resources/genes"

declare -A GTF=(
  [hg38]="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/genes.gtf"
  [hg19]="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/genes.gtf"
  [mm10]="/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/genes.gtf"
  [hs1]="/data/CCBR_Pipeliner/db/PipeDB/Indices/hs1/genes.gtf"
)

for g in hg38 hg19 mm10 hs1; do
  gtf="${GTF[$g]}"
  out="$OUT_ROOT/${g}.geneinfo.all.bed"
  out_pc="$OUT_ROOT/${g}.geneinfo.protein_coding.bed"

  echo ">> $g"
  awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
      gid=""; gname=""; gtype="";
      if (match($9,/gene_id "([^"]+)"/,m)) gid=m[1];
      if (match($9,/gene_name "([^"]+)"/,m)) gname=m[1]; else gname=gid;
      if (match($9,/(gene_type|gene_biotype) "([^"]+)"/,m)) gtype=m[2]; else gtype="NA";
      print $1, $4-1, $5, $7, gname, gtype, gid
    }' "$gtf" | sort -k1,1 -k2,2n > "$out"

  grep -w "protein_coding" "$out" > "$out_pc"
done
