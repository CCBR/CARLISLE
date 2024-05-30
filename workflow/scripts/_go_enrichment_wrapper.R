#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("--rmd", type="character", required=TRUE,
                    help="path to rmd")
parser$add_argument("--carlisle_functions", type="character", required=TRUE, 
          help="path to carlisle functions file")
parser$add_argument("--output_dir", type="character", required=FALSE,
                    help = "output_dir")
parser$add_argument("--report", type="character", required=TRUE,
                    help = "HTML report")
parser$add_argument("--peak_list", type="character", required=TRUE,
                    help = "peak_list")
parser$add_argument("--species", type="character", required=TRUE,
                    help = "species")
parser$add_argument("--geneset_id", type="character", required=FALSE, default="GOBP",
                    help = "geneset_id")
parser$add_argument("--dedup_status", type="character", required=FALSE, default="dedup",
                    help = "deduplication status")
args <- parser$parse_args()

debug="FALSE"
if (debug){
  carlisle_functions="/data/CCBR_Pipeliner/Pipelines/CARLISLE/latest/workflow/scripts/_carlisle_functions.R"
  output_dir="~/../../Volumes/data/tmp"
  report="~/../../Volumes/data/tmp/carlisle/report.html"
  peak_list="macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.narrow.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.narrow.peaks.bed"
  species="hg38"
  geneset_id="GOBP"
  dedup_status="dedup"
} else {
  carlisle_functions=args$carlisle_functions
  output_dir=args$output_dir
  report=args$report
  peak_list=args$peak_list
  species=args$species
  geneset_id=args$geneset_id
  dedup_status=args$dedup_status
}

parameters=list(
  carlisle_functions=carlisle_functions,
  output_dir=output_dir,
  peak_list=peak_list,
  species=species,
  geneset_id=geneset_id,
  dedup_status=dedup_status)

rmarkdown::render(args$rmd,
                  params=parameters,
                  output_file = report)