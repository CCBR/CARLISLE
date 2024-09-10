#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("--rmd",
  type = "character", required = TRUE,
  help = "path to rmd"
)
parser$add_argument("--carlisle_functions",
  type = "character", required = TRUE,
  help = "path to carlisle functions file"
)
parser$add_argument("--report",
  type = "character", required = TRUE,
  help = "HTML report"
)
parser$add_argument("--bam_list",
  type = "character", required = TRUE,
  help = "list of .idxstats bam files"
)
parser$add_argument("--spikein_control",
  type = "character", required = TRUE,
  help = "spike in species type"
)
args <- parser$parse_args()

debug <- "FALSE"
if (debug) {
  carlisle_functions <- "/data/CCBR_Pipeliner/Pipelines/CARLISLE/latest/workflow/scripts/_carlisle_functions.R"
  report <- "~/../../Volumes/data/tmp/carlisle/spike_report.html"
  bam_list <- "macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.narrow.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.narrow.peaks.bed"
  spikein_control <- "NC_000913.3"
} else {
  carlisle_functions <- args$carlisle_functions
  report <- args$report
  bam_list <- args$bam_list
  spikein_control <- args$spikein_control
}

parameters <- list(
  carlisle_functions = carlisle_functions,
  bam_list = bam_list,
  spikein_control = spikein_control
)

rmarkdown::render(args$rmd,
  params = parameters,
  output_file = report
)
