# Getting Started

## Introduction

The **[CARLISLE Pipeline](https://github.com/CCBR/CARLISLE)** is designed for the analysis of **CUT&RUN** and **CUT&Tag** datasets, two closely related chromatin profiling methods that provide high-resolution maps of protein-DNA interactions. Both techniques use antibody-targeted micrococcal nuclease (MNase) or Tn5 transposase to selectively cleave or tag DNA fragments bound by proteins of interest. The resulting DNA is then purified and sequenced to identify binding sites across the genome. CUT&RUN, in particular, offers improved signal-to-noise ratio and lower background compared to traditional **[ChIP-seq](https://www.nature.com/articles/nmeth.3327)**, making it highly efficient even with limited input material.

While CUT&Tag uses an engineered transposase to insert sequencing adapters directly into chromatin fragments, CUT&RUN relies on nuclease digestion and separate library preparation. Despite these differences, the downstream computational workflow—including alignment, filtering, deduplication, and peak calling—is largely similar. The CARLISLE pipeline supports both methods, offering flexibility for diverse chromatin profiling experiments and maintaining compatibility with established Biowulf HPC environments.

Inspired by the **[nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.2/)** pipeline, CARLISLE incorporates modular and reproducible analysis steps that emphasize transparency and scalability. It begins with raw **FASTQ** files, performing adapter trimming followed by alignment using **[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**. Linear deduplication is an important part of the process—especially for CUT&RUN and CUT&Tag data—to remove PCR artifacts and ensure accurate quantification of unique DNA fragments.

Normalization is performed using either user-provided spike-in controls (e.g., _E. coli_) or library size-based scaling. Peak calling is then executed using **[MACS2](https://github.com/macs3-project/MACS)**, **[SEACR](https://seacr.fredhutch.org/)**, and **[GoPeaks](https://github.com/maxsonBraunLab/gopeaks)**—with **GoPeaks** recommended for its robust performance on CUT&RUN data. Following peak calling, CARLISLE annotates results and summarizes them into detailed reports, with optional differential analysis handled via **[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)**. Quality control metrics are generated using **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** and **[MultiQC](https://multiqc.info/)**, while enrichment and annotation are supported through **[HOMER](http://homer.ucsd.edu/homer/)**, **[ROSE](https://bitbucket.org/young_computation/rose/src/master/)**, and **[ChIP-Enrich](https://chipenrich.med.umich.edu/)**.

---

## Setup Dependencies

> **Note:** CARLISLE supports execution via **[Singularity/Apptainer](https://apptainer.org/)** containers on **[Biowulf HPC](https://hpc.nih.gov/)**. Most dependencies are loaded from environment modules; select R-based steps use container images. Provide a local Singularity cache directory at run time with the `--singcache` flag.

CARLISLE relies on several dependencies, most of which are pre-installed and auto-loaded on **Biowulf** via the `ccbrpipeliner` module.

| Tool         | Module / Version                                          |
| ------------ | --------------------------------------------------------- |
| bedtools     | `bedtools/2.30.0`                                         |
| bedops       | `bedops/2.4.41`                                           |
| bowtie2      | `bowtie/2-2.4.5`                                          |
| cutadapt     | `cutadapt/4.0`                                            |
| deeptools    | `deeptools/3.5.1`                                         |
| fastqc       | `fastqc/0.11.9`                                           |
| fastq_screen | `fastq_screen/0.15.2`                                     |
| fastq_val    | `/data/CCBR_Pipeliner/bin/fastQValidator`                 |
| fastxtoolkit | `fastxtoolkit/0.0.14`                                     |
| gopeaks      | `/data/CCBR_Pipeliner/bin/GoPeaks/gopeaks`                |
| homer        | `homer/4.11.1`                                            |
| macs2        | `macs/2.2.7.1`                                            |
| meme         | `meme/5.5.5`                                              |
| multiqc      | `multiqc/1.14`                                            |
| perl         | `perl/5.34`                                               |
| picard       | `picard/2.27.3`                                           |
| python       | `python/3.9`                                              |
| R            | `R/4.3.2`                                                 |
| rose         | `ROSE/1.3.1`                                              |
| samtools     | `samtools/1.15`                                           |
| seacr        | `SEACR/1.4-beta.2`                                        |
| ucsc         | `ucsc/445`                                                |

---

## Login to the Cluster

CARLISLE has been tested and validated on **[NIH Biowulf HPC](https://hpc.nih.gov/)**.

```bash
# SSH into the Biowulf head node
ssh -Y $USER@biowulf.nih.gov
```

Move to your project directory before proceeding.

---

## Load an Interactive Session

An interactive session is required before executing any CARLISLE sub-commands, even if the job will later be submitted to the cluster.

```bash
# Request an interactive node
sinteractive --time=12:00:00 --mem=8gb --cpus-per-task=4 --pty bash
```

---

## Load Dependencies

Load the CARLISLE environment by activating the **ccbrpipeliner** module:

```bash
module load ccbrpipeliner/8
```

Check which version is currently active:

```bash
carlisle --version
[+] Loading singularity  4.2.2  on cn0001
[+] Loading snakemake  7.32.4
Pipeline Dir: /vf/users/CCBR_Pipeliner/Pipelines/CARLISLE/.v2.7.6
Version: 2.7.6
```

```
carlisle --help
[+] Loading singularity  4.2.2  on cn0001
[+] Loading snakemake  7.32.4
Pipeline Dir: /vf/users/CCBR_Pipeliner/Pipelines/CARLISLE/.v2.7.6
/spin1/home/linux/kopardevn/carlisle
  --> run CARLISLE
  Cut And Run anaLysIS pipeLinE

  USAGE:
    bash /vf/users/CCBR_Pipeliner/Pipelines/CARLISLE/.v2.7.6/carlisle -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
  Required Arguments:
  1.  RUNMODE: [Type: String] Valid options:
      *) init : initialize workdir
      *) run : run with slurm
      *) reset : DELETE workdir dir and re-init it
      *) dryrun : dry run snakemake to generate DAG
      *) unlock : unlock workdir if locked by snakemake
      *) runlocal : run without submitting to sbatch
      *) runtest: run on cluster with included hg38 test dataset
  2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.

  Optional Arguments:
      *) -f / --force: Force flag will re-initialize a previously initialized workdir
```
