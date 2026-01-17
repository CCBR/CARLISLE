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
parser$add_argument("--output_dir",
  type = "character", required = FALSE,
  help = "output_dir"
)
parser$add_argument("--report",
  type = "character", required = TRUE,
  help = "HTML report"
)
parser$add_argument("--contrast_file",
  type = "character", required = FALSE,
  help = "TSV produced by create_contrast_peakcaller_files (contains replicate IDs + peak BED paths)"
)
parser$add_argument("--peak_list",
  type = "character", required = FALSE,
  help = "DEPRECATED: xxx-delimited list of peak BED paths (kept for backward compatibility)"
)
parser$add_argument("--species",
  type = "character", required = TRUE,
  help = "species"
)
parser$add_argument("--geneset_id",
  type = "character", required = FALSE, default = "GOBP",
  help = "geneset_id"
)
parser$add_argument("--dedup_status",
  type = "character", required = FALSE, default = "dedup",
  help = "deduplication status"
)
parser$add_argument("--n_cores",
  type = "integer", required = FALSE, default = 1,
  help = "Number of CPU cores to use inside chipenrich/polyonrich/broadenrich/hybridenrich where supported"
)
parser$add_argument("--skip_hybrid",
  action = "store_true",
  help = "Skip the hybrid enrichment (peak_h / hybridenrich) step. Useful for quick testing."
)
args <- parser$parse_args()

debug <- "FALSE"
if (debug) {
  carlisle_functions <- "/data/CCBR_Pipeliner/Pipelines/CARLISLE/latest/workflow/scripts/_carlisle_functions.R"
  output_dir <- "~/../../Volumes/data/tmp"
  report <- "~/../../Volumes/data/tmp/carlisle/report.html"
  peak_list <- "macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_1_vs_nocontrol.dedup.narrow.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.broad.peaks.bed macs2/peak_output/53_H3K4me3_2_vs_nocontrol.dedup.narrow.peaks.bed"
  contrast_file <- NULL
  species <- "hg38"
  geneset_id <- "GOBP"
  dedup_status <- "dedup"
  n_cores <- 1
  skip_hybrid <- FALSE
} else {
  carlisle_functions <- args$carlisle_functions
  output_dir <- args$output_dir
  report <- args$report
  peak_list <- args$peak_list
  contrast_file <- args$contrast_file
  species <- args$species
  geneset_id <- args$geneset_id
  dedup_status <- args$dedup_status
  n_cores <- args$n_cores
  skip_hybrid <- isTRUE(args$skip_hybrid)
}

parameters <- list(
  carlisle_functions = carlisle_functions,
  output_dir = output_dir,
  peak_list = peak_list,
  contrast_file = contrast_file,
  species = species,
  geneset_id = geneset_id,
  dedup_status = dedup_status,
  n_cores = n_cores,
  skip_hybrid = skip_hybrid
)

rmarkdown::render(args$rmd,
  params = parameters,
  output_file = report
)
