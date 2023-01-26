`gtf2tssBed.py` script is used to generate the BED file from GTF annotations.
```bash
./gtf2tssBed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/genes.gtf --bed hg38.tss.bed
./gtf2tssBed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/genes.gtf --bed hg19.tss.bed
./gtf2tssBed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/T2T/genes.gtf --bed T2T.tss.bed
```
> NOTE: Previously, the following files were used as the tss BEDs but these:
> * are 1-based (BEDs are 0-based)
> * have 5 columns (BEDs have 3 or 6 or 12 columns)
> * column 5 is strand (BEDs have column 6 reserved for strand)
> 
>  Hence, we have replaced these with the new BED files created using the above commands.
> #### hg38.tss.bed source
> https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/ref/hg38.tss.bed
> #### hg19.tss.bed source
> https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/ref/hg19.tss.bed
