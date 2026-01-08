#!/usr/bin/env bash

# make sure that all files are decompressed before running this!

set -euo pipefail

cd /data/CCBR_Pipeliner/Pipelines/CARLISLE/vnk_dev/resources/cCREs

for g in hg38 hs1 hg19 mm10; do
  infile="${g}-cCREs.bed"
  if [[ ! -f "$infile" ]]; then
    echo "Missing $infile; skipping"
    continue
  fi

  echo ">> Splitting $infile"
  rm -f ${g}.*.bed
  awk -v OFS='\t' -v prefix="${g}." '{
    f = prefix $NF ".bed";
    print >> f;
    files[f]=1
  }
  END { for (f in files) close(f) }' "$infile"
done
