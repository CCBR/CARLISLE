#!/bin/bash

module load ucsc
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons /data/CCBR_Pipeliner/db/PipeDB/Indices/hs1/genes.gtf hs1.genes.genePred
head -n1 hg19_refseq.ucsc > hs1.ucsc
awk -F"\t" -v OFS="\t" '{print $NF,$_}' hs1.genes.genePred >> hs1.ucsc