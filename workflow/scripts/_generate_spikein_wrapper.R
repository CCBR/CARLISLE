#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("--rmd", type="character", required=TRUE,
                    help="path to rmd")
parser$add_argument("--sourcefile", type="character", required=TRUE,
                    help="path to function file")
parser$add_argument("--tmpdir", type="character", required=FALSE, default="/tmp",
                    help = "tmpdir")
parser$add_argument("--report", type="character", required=TRUE,
                    help = "HTML report")
parser$add_argument("--bam_list", type="character", required=TRUE,
                    help = "list of .idxstats bam files")
parser$add_argument("--spikein_control", type="character", required=TRUE,
                    help = "spike in species type")
args <- parser$parse_args()

debug="FALSE"
if (debug){
  sourcefile="~/../../Volumes/Pipelines/CARLISLE_dev/workflow/scripts/_carlisle_functions.R"
  tmpdir="/dev/shm"
  report="~/../../Volumes/data/tmp/carlisle/spike_report.html"
  bam_list="macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.narrow.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.narrow.peaks.bed"
  spikein_control="NC_000913.3"
} else {
  sourcefile=args$sourcefile
  tmpdir=args$tmpdir
  report=args$report
  bam_list=args$bam_list
  spikein_control=args$spikein_control
}

parameters=list(sourcefile=sourcefile,
                bam_list=bam_list,
                spikein_control=spikein_control)

rmarkdown::render(args$rmd,
                  params=parameters,
                  output_file = report,
                  intermediates_dir = paste(tmpdir,"intermediates_dir",sep="/"),
                  knit_root_dir = paste(tmpdir,"knit_root_dir",sep="/"))