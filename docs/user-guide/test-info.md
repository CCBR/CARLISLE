# Pipeline Tutorial

Welcome to the **CARLISLE Pipeline Tutorial!** This guide walks you through running the CARLISLE pipeline using the provided test dataset on the **[NIH Biowulf HPC](https://hpc.nih.gov/)** environment.

---

## Getting Started

Before beginning, review the [Getting Started Guide](https://ccbr.github.io/CARLISLE/user-guide/getting-started/) for installation, environment setup, and dependency loading instructions.


---

### Step 1. Set Your Working Directory

Navigate to your project directory on Biowulf:

```bash
cd /path/to/your/project/directory
```

---

### Step 2. Initialize the Pipeline

Load the CARLISLE module and initialize your working directory:

```bash
module load ccbrpipeliner
carlisle --runmode=init --workdir=/path/to/output/dir
```

This command copies the required configuration, manifest, and Snakefiles into your chosen output directory (`WORKDIR`). Initialization must be done before any other CARLISLE operation.

---

## Submitting the Test Data

The test dataset provided with CARLISLE enables you to validate the installation and confirm correct execution. The test includes minimal FASTQ files, configurations, and manifests.

### Step 3. Run the Test Command

Execute the built-in test run to validate pipeline functionality:

```bash
carlisle --runmode=runtest --workdir=/path/to/output/dir
```

This command prepares the test data, performs a dry-run to validate workflow dependencies, and then submits the pipeline to the Biowulf SLURM cluster.

---

### Expected Output

During a successful test run, you should see a job summary similar to the one below, detailing the number of tasks executed per Snakemake rule:

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

> 💡 **Tip:** This job summary confirms successful rule execution, resource allocation, and workflow orchestratio
