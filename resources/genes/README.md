## GTF

- hg38 = /data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/genes.gtf
- hg19 = /data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/genes.gtf
- mm10 = /data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/mm10/gencode.vM25.annotation.gtf
- hs1 = /data/CCBR_Pipeliner/db/PipeDB/Indices/hs1/genes.gtf

## sizes

- hg38: `cut -f1,2 /data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/hg38.fa.fai > hg38.sizes`
- hg19: `cut -f1,2 /data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/ref.fa.fai > hg19.sizes`
- mm10: `cut -f1,2 /data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/ref.fa.fai > mm10.sizes`
- hs1:  `cp /data/CCBR_Pipeliner/db/PipeDB/Indices/hs1/hs1.sizes T2T.sizes`

## ambiguous BED

`bash get_ambiguous_bed.sh``

