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
parser$add_argument("--countsmatrix", type="character", required=TRUE,
                    help="countsmatrix as TSV")
parser$add_argument("--sampleinfo", type="character", required=TRUE,
                    help="sample info as TSV")
parser$add_argument("--dupstatus", type="character", required=TRUE,
                    help="either dedup or no_dedup")
parser$add_argument("--cpm", type="double", default=1, required=FALSE,
                    help="cpm cutoff [default %(default)s]")
parser$add_argument("--condition1", type="character", required=TRUE,
                    help = "condition1")
parser$add_argument("--condition2", type="character", required=TRUE,
                    help = "condition2")
parser$add_argument("--results", type="character", required=TRUE,
                    help = "path to results TSV")
parser$add_argument("--report", type="character", required=TRUE,
                    help = "HTML report")


args <- parser$parse_args()

if (debug){
  rawcountsmatrix="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/counts_matrix.txt"
  coldata="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/sample_info.txt"
  dupstatus="dedup"
  condition1="siSmyd3_2m_Smyd3_0.25HCHO_500K"
  condition2="siNC_2m_Smyd3_0.25HCHO_500K"
  indexcols="peakID"
  cpm_cutoff=1
  results="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/results.txt"
  report="~/CCBR/projects/ccbr1155/CS030586_CARAP/diff/report.html"
} else {
  rawcountsmatrix=args$countsmatrix
  coldata=args$sampleinfo
  dupstatus=args$dupstatus
  condition1=args$condition1
  condition2=args$condition2
  cpm_cutoff=args$cpm
  indexcols="peakID"
  results=args$results
  report=args$report
}


parameters=list(rawcountsmatrix=rawcountsmatrix,
            coldata=coldata,
            dupstatus=dupstatus,
            condition1=condition1,
            condition2=condition2,
            indexcols=indexcols,
            cpm_cutoff=cpm_cutoff,
            results=results)

rmarkdown::render(args$rmd,
  params=parameters,
  output_file = report)