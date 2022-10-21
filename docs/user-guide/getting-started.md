# Overview
The CARLISLE github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

## 1. Getting Started
1.1 Introduction
The CARLISLE Pipelie beings with raw FASTQ files and performs trimming followed by alignment using [BOWTIE2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Data is then normalized through either the use of an user-species species (IE E.Coli) spike-in control or through the determined library size. Peaks are then called using [SEACR](https://github.com/FredHutch/SEACR) and [MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html), with various options selected by the user. Peaks are then annotated, and summarized into reports. If designated, differential analysis is performed using [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). QC reports are also generated with each project.

The following are sub-commands used within CARLISLE:

- initialize: initalize the pipeline
- dryrun: predict the binding of peptides to any MHC molecule
- cluster: execute the pipeline on the Biowulf HPC
- local: execute a local, interactive, session
- git: execute GitHub actions
- unlock: unlock directory
- DAG: create DAG report
- report: create SNAKEMAKE report

## 1.2 Setup Dependencies
CARLISLE has several dependencies listed below. These dependencies can be installed by a sysadmin. All dependencies will be automatically loaded if running from Biowulf.

- bedtools: "bedtools/2.30.0"
- bedops: "bedops/2.4.40"
- bowtie2: "bowtie/2-2.4.2"
- cutadapt: "cutadapt/1.18"
- fastqc: "fastqc/0.11.9"
- fastxtoolkit: "fastxtoolkit/0.0.14"
- macs2: "macs/2.2.7.1"
- multiqc: "multiqc/1.9"
- picard: "picard/2.26.9"
- python37: "python/3.7"
- R: "R/4.0"
- samtools: "samtools/1.15"
- seacr: "seacr/1.4-beta.2"
- ucsc: "ucsc/407"

## 1.3 Login to the cluster
CARLISLE has been exclusively tested on Biowulf HPC. Login to the cluster's head node and move into the pipeline location.
```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov
```

## 1.4 Load an interactive session 
An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.

```
# Grab an interactive node
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
```
