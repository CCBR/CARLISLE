# Overview

The CARLISLE github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

## 1. Getting Started

## 1.1 Introduction

The CARLISLE Pipelie beings with raw FASTQ files and performs trimming followed by alignment using [BOWTIE2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Data is then normalized through either the use of an user-species species (IE E.Coli) spike-in control or through the determined library size. Peaks are then called using [MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html), [SEACR](https://github.com/FredHutch/SEACR), and [GoPEAKS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02707-w) with various options selected by the user. Peaks are then annotated, and summarized into reports. If designated, differential analysis is performed using [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). QC reports are also generated with each project using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MULTIQC](https://multiqc.info/). Annotations are added using [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html) and [ROSE](https://github.com/stjude/ROSE). GSEA Enrichment analysis predictions are added using [CHIPENRICH](https://bioconductor.org/packages/devel/bioc/vignettes/chipenrich/inst/doc/chipenrich-vignette.html).

The following are sub-commands used within CARLISLE:

- initialize: initalize the pipeline
- dryrun: predict the binding of peptides to any MHC molecule
- cluster: execute the pipeline on the Biowulf HPC
- local: execute a local, interactive, session
- git: execute GitHub actions
- unlock: unlock directory
- DAG: create DAG report
- report: create SNAKEMAKE report
- runtest: copies test manifests and files to WORKDIR

## 1.2 Setup Dependencies

CARLISLE has several dependencies listed below. These dependencies can be installed by a sysadmin. All dependencies will be automatically loaded if running from Biowulf.

- bedtools: "bedtools/2.30.0"
- bedops: "bedops/2.4.40"
- bowtie2: "bowtie/2-2.4.2"
- cutadapt: "cutadapt/1.18"
- fastqc: "fastqc/0.11.9"
- fastq_screen: "fastq_screen/0.15.2"
- fastq_val: "/data/CCBR_Pipeliner/iCLIP/bin/fastQValidator"
- fastxtoolkit: "fastxtoolkit/0.0.14"
- gopeaks: "github clone https://github.com/maxsonBraunLab/gopeaks"
- macs2: "macs/2.2.7.1"
- multiqc: "multiqc/1.9"
- perl: "perl/5.34.0"
- picard: "picard/2.26.9"
- python37: "python/3.7"
- R: "R/4.2.2"
- rose: "ROSE/1.3.1"
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
sinteractive --time=12:00:00 --mem=8gb  --cpus-per-task=4 --pty bash
```
