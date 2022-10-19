# 4. Expected Outputs
The following directories are created under the WORKDIR/results directory:

- alignment_stats: this directory include information on the alignment of each sample
- bigwig: this directory includes the bigwig files for each sample 
- peaks:
    - contrasts: this includes the contrasts for each line listed in the contrast manifest
    - seacr: this includes all peak calls from SEACR for each sample
    - macs2: this includes all peak calls from MACS2 for each sample
- bam: this includes BAM files, statistics on samples, statistics on spike-in controls for each sample
- bedgraph: this includes BEDGRAPH files and statistic summaries for each sample