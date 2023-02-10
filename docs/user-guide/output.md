# 4. Expected Outputs
The following directories are created under the WORKDIR/results directory:

- annotation: this directory will include annotations from HOMER and ROSE. It contains a sub-directory that relates to the quality threshold used.
- alignment_stats: this directory include information on the alignment of each sample
- bigwig: this directory includes the bigwig files for each sample 
- peaks: It contains a sub-directory that relates to the quality threshold used.
    - contrasts: this includes the contrasts for each line listed in the contrast manifest
    - seacr: this includes all peak calls from SEACR for each sample
    - macs2: this includes all peak calls from MACS2 for each sample
    - gopeaks: this includes all peak calls from GOPEAKS for each sample
- bam: this includes BAM files, statistics on samples, statistics on spike-in controls for each sample
- bedgraph: this includes BEDGRAPH files and statistic summaries for each sample