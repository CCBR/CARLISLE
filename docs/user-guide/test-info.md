# 5. Pipeline Tutorial
Welcome to the CARLISLE Pipeline Tutorial!

## 5.1 Getting Started
Review the information on the [Getting Started](https://ccbr.github.io/CARLISLE/user-guide/getting-started/) for a complete overview the pipeline. The tutorial below will use test data available on NIH Biowulf HPC only. All example code will assume you are running v1.0 of the pipeline, using test data available on GitHub.

A. Change working directory to the CARLISLE repository

B. Initialize Pipeline
```
bash ./path/to/dir/carlisle --runmode=init --workdir=/path/to/output/dir
```

## 5.2 Submit the test data
Test data is included in the .test directory as well as the config directory.

A Run the test command to prepare the data, perform a dry-run and submit to the cluster
```
bash ./path/to/dir/carlisle --runmode=testrun --workdir=/path/to/output/dir

```

- An expected output for the `testrun` is as follows:
```
Job stats:
job                              count    min threads    max threads
-----------------------------  -------  -------------  -------------
DESeq                                  24              1              1
align                                   9             56             56
alignstats                              9              2              2
all                                     1              1              1
bam2bg                                  9              4              4
create_contrast_data_files             24              1              1
create_contrast_peakcaller_files       12              1              1
create_reference                        1             32             32
create_replicate_sample_table           1              1              1
diffbb                                 24              1              1
filter                                 18              2              2
findMotif                              96              6              6
gather_alignstats                       1              1              1
go_enrichment                          12              1              1
gopeaks_broad                          16              2              2
gopeaks_narrow                         16              2              2
macs2_broad                            16              2              2
macs2_narrow                           16              2              2
make_counts_matrix                     24              1              1
multiqc                                 2              1              1
qc_fastqc                               9              1              1
rose                                   96              2              2
seacr_relaxed                          16              2              2
seacr_stringent                        16              2              2
spikein_assessment                      1              1              1
trim                                    9             56             56
total                                 478              1             56
```

## 5.3 Review outputs
Review the expected outputs on the [Output](https://ccbr.github.io/CARLISLE/user-guide/output/) page. If there are errors, review and performing stesp described on the [Troubleshooting](https://ccbr.github.io/CARLISLE/user-guide/troubleshooting/) page as needed.