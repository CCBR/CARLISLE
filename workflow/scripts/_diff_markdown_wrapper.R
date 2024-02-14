#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

#debug="TRUE"
#path="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff"

debug="FALSE"
path="/data/CCBR/projects/ccbr1155/CS030586_CARAP/diff"

# create parser object
parser <- ArgumentParser()
parser$add_argument("--rmd", type="character", required=TRUE,
          help="path to rmd")
parser$add_argument("--carlisle_functions", type="character", required=TRUE, 
          help="path to carlisle functions file")
parser$add_argument("--Rlib_dir", type="character", required=TRUE, 
          help="path to R lib directory")
parser$add_argument("--Rpkg_config", type="character", required=TRUE, 
          help="path to package config")
parser$add_argument("--spiked", type="character", required=TRUE, 
          help="type of normalization used")
parser$add_argument("--rawcountsprescaled", action='store_true',
          help="if counts are scaled by spike-in already ... Y (for AUC-based method) or N (for fragments-based method)")
parser$add_argument("--scalesfbymean", action='store_true',
          help="DESeq2 scaling factors are around 1. To ensure that spike-in scaling factors are also around 1 divide each scaling factor by mean of all scaling factors.")
parser$add_argument("--htsfilter", action='store_true',
          help="Use HTSFilter")
parser$add_argument("--contrast_data", type="character", required=FALSE, default=NULL,
          help="contrast_data inputs file")
parser$add_argument("--countsmatrix", type="character", required=TRUE,
          help="countsmatrix as TSV")
parser$add_argument("--sampleinfo", type="character", required=TRUE,
          help="sample info as TSV")
parser$add_argument("--dupstatus", type="character", required=TRUE,
          help="either dedup or no_dedup")
parser$add_argument("--fdr_cutoff", type="double", default=0.05, required=FALSE,
          help="FDR cutoff [default %(default)s]")
parser$add_argument("--log2fc_cutoff", type="double", default=0.59, required=FALSE,
          help="log2foldchange cutoff [default %(default)s]")
parser$add_argument("--condition1", type="character", required=TRUE,
          help = "condition1")
parser$add_argument("--condition2", type="character", required=TRUE,
          help = "condition2")
parser$add_argument("--results", type="character", required=TRUE,
          help = "path to results TSV")
parser$add_argument("--report", type="character", required=TRUE,
          help = "HTML report")
parser$add_argument("--elbowlimits", type="character", required=TRUE,
          help = "YAML ELBOW limits")
parser$add_argument("--tmpdir", type="character", required=FALSE, default="/tmp",
          help = "tmpdir")
parser$add_argument("--species", type="character", required=TRUE,
          help = "species")
parser$add_argument("--gtf", type="character", required=FALSE,
          help = "gtf path - needed for HS1")
args <- parser$parse_args()

gtf <- args$gtf
if (debug){
  carlisle_functions="/data/CCBR_Pipeliner/Pipelines/CARLISLE/latest/workflow/scripts/_carlisle_functions.R"
  Rlib_dir="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_4.3_carlisle/"
  Rpkg_config="/data/CCBR_Pipeliner/Pipelines/CARLISLE/latest/conf/Rpack.config"
  rawcountsmatrix="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/counts_matrix.txt"
  coldata="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/sample_info.txt"
  dupstatus="dedup"
  condition1="siSmyd3_2m_Smyd3_0.25HCHO_500K"
  condition2="siNC_2m_Smyd3_0.25HCHO_500K"
  indexcols="peakID"
  fdr_cutoff=0.05
  log2fc_cutoff=0.59
  results="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/results.txt"
  report="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/report.html"
  spiked="LIBRARY"
  rawcountsprescaled="N"
  scalesfbymean="N"
  htsfilter="Y"
  elbowlimits="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/elbow.yaml"
  tmpdir="/dev/shm"
  gtf="~/../../Volumes/CCBR_Pipeliner/db/PipeDB/Indices/hs1/genes.gtf"
} else {
  carlisle_functions=args$carlisle_functions
  Rlib_dir=args$Rlib_dir
  Rpkg_config=args$Rpkg_config
  rawcountsmatrix=args$countsmatrix
  coldata=args$sampleinfo
  dupstatus=args$dupstatus
  condition1=args$condition1
  condition2=args$condition2
  fdr_cutoff=args$fdr_cutoff
  log2fc_cutoff=args$fdr_cutoff
  indexcols="peakID"
  results=args$results
  report=args$report
  spiked=args$spiked
  contrast_data=args$contrast_data
  elbowlimits=args$elbowlimits
  if (args$rawcountsprescaled) {rawcountsprescaled="Y"} else {rawcountsprescaled="N"}
  if (args$scalesfbymean) {scalesfbymean="Y"} else {scalesfbymean="N"}
  tmpdir=args$tmpdir
  if (args$htsfilter) {htsfilter="Y"} else {htsfilter="N"}
  species=args$species
  gtf=args$gtf
}

parameters=list(
  carlisle_functions=carlisle_functions,
  Rlib_dir=Rlib_dir,
  Rpkg_config=Rpkg_config,
  rawcountsmatrix=rawcountsmatrix,
  coldata=coldata,
  spiked=spiked,
  rawcountsprescaled=rawcountsprescaled,
  scalesfbymean=scalesfbymean,
  contrast_data=contrast_data,
  dupstatus=dupstatus,
  condition1=condition1,
  condition2=condition2,
  indexcols=indexcols,
  htsfilter=htsfilter,
  fdr_cutoff=fdr_cutoff,
  log2fc_cutoff=log2fc_cutoff,
  results=results,
  elbowlimits=elbowlimits,
  species=species,
  gtf=gtf)

rmarkdown::render(args$rmd,
  params=parameters,
  output_file = report,
  intermediates_dir = paste(tmpdir,"intermediates_dir",sep="/"),
  knit_root_dir = paste(tmpdir,"knit_root_dir",sep="/"))