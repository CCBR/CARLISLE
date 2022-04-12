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
parser$add_argument("--spiked", type="character", required=TRUE, 
                    help="if spike-in is present or not ... Y or N")
parser$add_argument("--rawcountsprescaled", action='store_true',
                    help="if counts are scaled by spike-in already ... Y (for AUC-based method) or N (for fragments-based method)")
parser$add_argument("--scalesfbymean", action='store_true',
                    help="DESeq2 scaling factors are around 1. To ensure that spike-in scaling factors are also around 1 divide each scaling factor by mean of all scaling factors.")
parser$add_argument("--htsfilter", action='store_true',
                    help="Use HTSFilter")
parser$add_argument("--bbpaths", type="character", required=FALSE, default=NULL,
                    help="bedbedgraph inputs file")
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


args <- parser$parse_args()

if (debug){
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
  spiked="N"
  rawcountsprescaled="N"
  scalesfbymean="N"
  htsfilter="Y"
  elbowlimits="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/elbow.yaml"
  tmpdir="/dev/shm"
} else {
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
  bbpaths=args$bbpaths
  elbowlimits=args$elbowlimits
  if (args$rawcountsprescaled) {rawcountsprescaled="Y"} else {rawcountsprescaled="N"}
  if (args$scalesfbymean) {scalesfbymean="Y"} else {scalesfbymean="N"}
  tmpdir=args$tmpdir
  if (args$htsfilter) {htsfilter="Y"} else {htsfilter="N"}
}


parameters=list(rawcountsmatrix=rawcountsmatrix,
            coldata=coldata,
            spiked=spiked,
            rawcountsprescaled=rawcountsprescaled,
            scalesfbymean=scalesfbymean,
            bbpaths=bbpaths,
            dupstatus=dupstatus,
            condition1=condition1,
            condition2=condition2,
            indexcols=indexcols,
            htsfilter=htsfilter,
            fdr_cutoff=fdr_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            results=results,
            elbowlimits=elbowlimits)

rmarkdown::render(args$rmd,
  params=parameters,
  output_file = report,
  intermediates_dir = paste(tmpdir,"intermediates_dir",sep="/"),
  knit_root_dir = paste(tmpdir,"knit_root_dir",sep="/"))