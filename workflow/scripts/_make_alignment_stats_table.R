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
parser$add_argument("--scaleConstant",
  type = "double", default = 1e6, required = FALSE,
  help = "scaling constant, typically 1e6.[default %(default)s]"
)
parser$add_argument("--outTable",
  type = "character", required = TRUE,
  help = "absolute path to output table TSV file."
)

args <- parser$parse_args()
debug <- 0
if (debug == 1) {
  yaml_dir <- "/home/kopardevn/CCBR/projects/ccbr1155/CS030586_CARAP/results/alignment_stats"
  exclude_from_name <- ".alignment_stats.yaml"
  scale_constant <- 1e6
  out_table <- "test.tsv"
} else {
  yaml_dir <- args$yamlDir
  exclude_from_name <- args$excludeFromName
  scale_constant <- args$scaleConstant
  out_table <- args$outTable
}
setwd(yaml_dir)
files <- list.files(
  path = yaml_dir,
  pattern = "*.yaml"
)
df <- data.frame()
for (f in files) {
  sampleName <- gsub(pattern = exclude_from_name, replacement = "", f)
  readinyaml <- as.data.frame(read_yaml(f, fileEncoding = "UTF-8"))
  rownames(readinyaml) <- c(sampleName)
  df <- rbind(df, readinyaml)
}
df$scaling_factor <- scale_constant / df$dedup_nreads_spikein
df$duplication_rate_genome <- (df$raw_nreads_genome - df$dedup_nreads_genome) / df$raw_nreads_genome
df$duplication_rate_spikein <- (df$raw_nreads_spikein - df$dedup_nreads_spikein) / df$raw_nreads_spikein
column_order <- c(
  "sample_name",
  "nreads",
  "raw_nreads_genome",
  "raw_nreads_spikein",
  "scaling_factor",
  "duplication_rate_genome",
  "duplication_rate_spikein",
  "no_dedup_nreads_genome",
  "no_dedup_nreads_spikein",
  "dedup_nreads_genome",
  "dedup_nreads_spikein"
)


df %>% rownames_to_column(var = "sample_name") -> df
df <- df[, column_order]

write.table(df,
  file = out_table,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
