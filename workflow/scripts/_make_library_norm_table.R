#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("tidyverse"))


# create parser object
parser <- ArgumentParser()


parser$add_argument("--yamlDir",
  type = "character", required = TRUE,
  help = "absolute path to location of yamls"
)
parser$add_argument("--excludeFromName",
  type = "character", required = TRUE,
  help = "sample info as TSV"
)
parser$add_argument("--outTable",
  type = "character", required = TRUE,
  help = "absolute path to output table TSV file."
)
args <- parser$parse_args()

# assign params variables
debug <- 0
if (debug == 1) {
  yaml_dir <- "~/../../Volumes/CUTRUN/analysis/CS030666/carlisle_230118/results/alignment_stats"
  exclude_from_name <- ".alignment_stats.yaml"
  scale_constant <- 1e6
  out_table <- "test.tsv"
} else {
  yaml_dir <- args$yamlDir
  exclude_from_name <- args$excludeFromName
  scale_constant <- args$scaleConstant
  out_table <- args$outTable
}
# grab all files
setwd(yaml_dir)
files <- list.files(
  path = yaml_dir,
  pattern = "*.yaml"
)

# read in df and prep
df <- data.frame()
for (f in files) {
  sampleName <- gsub(pattern = exclude_from_name, replacement = "", f)
  readinyaml <- as.data.frame(read_yaml(f, fileEncoding = "UTF-8"))
  rownames(readinyaml) <- c(sampleName)
  df <- rbind(df, readinyaml)
}

# run once for dedup, once for nondedup
final_df <- data.frame()
for (dedup_type in c("dedup_nreads_genome", "no_dedup_nreads_genome")) {
  # determine scaling factor dependent on library size
  col_median <- median(df[, dedup_type])
  if (col_median > 100000000) {
    lib_factor <- 1e8
  } else if (col_median > 10000000) {
    lib_factor <- 1e7
  } else if (col_median > 1000000) {
    lib_factor <- 1e6
  } else if (col_median > 100000) {
    lib_factor <- 1e5
  } else if (col_median > 10000) {
    lib_factor <- 1e4
  } else if (col_median > 1000) {
    lib_factor <- 1e3
  } else if (col_median > 100) {
    lib_factor <- 1e2
  } else {
    lib_factor <- 1e1
  }

  df$library_size <- lib_factor / df[, dedup_type]
  df$sampleid <- rownames(df)
  df$dedup_type <- strsplit(dedup_type, "_")[[1]][1]
  df$dedup_type <- gsub("no", "no_dedup", df$dedup_type)

  # create final df
  select_cols <- c("sampleid", "library_size", "dedup_type")
  final_df <- rbind(final_df, df[, select_cols])
}

write.table(final_df, out_table, row.names = FALSE)
