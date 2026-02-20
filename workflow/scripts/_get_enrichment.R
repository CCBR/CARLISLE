#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--peaks_bed",
  type = "character", required = TRUE,
  help = "Path to peaks BED file"
)
parser$add_argument("--geneset_id",
  type = "character", required = TRUE,
  help = "GO geneset id(s), comma-separated (e.g., GOBP,GOCC,GOMF)"
)
parser$add_argument("--genome",
  type = "character", required = TRUE,
  help = "Genome: hg19, hg38, or mm10"
)
parser$add_argument("--output_tsv",
  type = "character", required = TRUE,
  help = "TSV output path for enrichment results"
)
parser$add_argument("--locusdef",
  type = "character", required = FALSE, default = "nearest_tss",
  help = "Locus definition for chipenrich"
)
parser$add_argument("--methods",
  type = "character", required = FALSE, default = "chipenrich",
  help = "Comma-separated methods: chipenrich,polyenrich,hybridenrich"
)
parser$add_argument("--n_cores",
  type = "integer", required = FALSE, default = 1,
  help = "Number of CPU cores to use where supported"
)
args <- parser$parse_args()

allowed_genomes <- c("hg19", "hg38", "mm10")
if (!(args$genome %in% allowed_genomes)) {
  stop(paste0("Invalid genome: ", args$genome, ". Must be one of: ", paste(allowed_genomes, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(chipenrich)
  library(AnnotationDbi)
})

output_dir <- dirname(args$output_tsv)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!is.null(args$n_cores)) {
  options(mc.cores = as.integer(args$n_cores))
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (isTRUE(args$n_cores > 1)) {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = as.integer(args$n_cores)), default = TRUE)
    } else {
      BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    }
  }
}

read_peak_file <- function(peak_file_in) {
  peak_df <- read.csv(peak_file_in, sep = "\t", header = FALSE)[, c("V1", "V2", "V3")]
  colnames(peak_df) <- c("chrom", "start", "end")
  peak_df
}

sanitize_tag <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) {
    x <- "geneset"
  }
  x
}

peaks_df <- read_peak_file(args$peaks_bed)

geneset_ids <- trimws(unlist(strsplit(args$geneset_id, ",")))
geneset_ids <- geneset_ids[nzchar(geneset_ids)]

methods <- tolower(trimws(unlist(strsplit(args$methods, ","))))
methods <- methods[nzchar(methods)]
allowed_methods <- c("chipenrich", "polyenrich", "hybridenrich")
if (!all(methods %in% allowed_methods)) {
  bad <- methods[!(methods %in% allowed_methods)]
  stop(paste0("Invalid method(s): ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed_methods, collapse = ", ")))
}
if (length(methods) == 0) {
  methods <- "chipenrich"
}

if (nrow(peaks_df) == 0) {
  empty_cols <- c(
    "Geneset.Type", "Geneset.ID", "Description", "P.value", "FDR", "Effect",
    "Odds.Ratio", "Status", "N.Geneset.Genes", "N.Geneset.Peak.Genes",
    "Geneset.Avg.Gene.Length", "Geneset.Peak.Genes"
  )
  empty_df <- as.data.frame(setNames(replicate(length(empty_cols), logical(0), simplify = FALSE), empty_cols))
  for (gs in geneset_ids) {
    gs_tag <- sanitize_tag(gs)
    for (m in methods) {
      out_tsv <- args$output_tsv
      out_tsv <- sub("\\.tsv$", paste0(".", gs_tag, ".", m, ".tsv"), out_tsv)
      if (!grepl("\\.tsv$", out_tsv)) {
        out_tsv <- paste0(out_tsv, ".", gs_tag, ".", m, ".tsv")
      }
      write.table(empty_df, file = out_tsv, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  }
  quit(status = 0)
}

for (gs in geneset_ids) {
  for (m in methods) {
    if (m == "chipenrich") {
      results <- chipenrich::chipenrich(
        peaks = peaks_df,
        genome = args$genome,
        genesets = gs,
        locusdef = args$locusdef,
        qc_plots = FALSE,
        out_name = NULL,
        n_cores = as.integer(args$n_cores)
      )
    } else if (m == "polyenrich") {
      results <- chipenrich::polyenrich(
        peaks = peaks_df,
        genome = args$genome,
        genesets = gs,
        locusdef = args$locusdef,
        method = "polyenrich",
        qc_plots = FALSE,
        out_name = NULL,
        n_cores = as.integer(args$n_cores)
      )
    } else {
      results <- chipenrich::hybridenrich(
        peaks = peaks_df,
        genome = args$genome,
        genesets = gs,
        locusdef = args$locusdef,
        qc_plots = FALSE,
        out_name = NULL,
        n_cores = as.integer(args$n_cores),
        min_geneset_size = 10
      )
    }

    result_out <- results$results

    if ("Geneset.Peak.Genes" %in% names(result_out)) {
      entrez_ids <- unique(unlist(strsplit(paste(result_out$Geneset.Peak.Genes, collapse = ","), ",")))
      entrez_ids <- trimws(entrez_ids)
      entrez_ids <- entrez_ids[nzchar(entrez_ids)]

      if (length(entrez_ids) > 0) {
        if (args$genome %in% c("hg19", "hg38")) {
          suppressPackageStartupMessages(library(org.Hs.eg.db))
          sym_map <- AnnotationDbi::select(org.Hs.eg.db,
            keys = entrez_ids, keytype = "ENTREZID", columns = c("SYMBOL")
          )
        } else {
          suppressPackageStartupMessages(library(org.Mm.eg.db))
          sym_map <- AnnotationDbi::select(org.Mm.eg.db,
            keys = entrez_ids, keytype = "ENTREZID", columns = c("SYMBOL")
          )
        }

        id_to_symbol <- sym_map %>%
          dplyr::filter(!is.na(SYMBOL)) %>%
          dplyr::distinct(ENTREZID, SYMBOL) %>%
          dplyr::group_by(ENTREZID) %>%
          dplyr::summarise(SYMBOL = SYMBOL[1], .groups = "drop")

        symbol_lookup <- setNames(id_to_symbol$SYMBOL, id_to_symbol$ENTREZID)

        convert_ids <- function(x) {
          ids <- trimws(unlist(strsplit(x, ",")))
          ids <- ids[nzchar(ids)]
          syms <- unname(symbol_lookup[ids])
          syms <- ifelse(is.na(syms) | !nzchar(syms), ids, syms)
          paste(syms, collapse = ", ")
        }

        result_out$Geneset.Peak.Genes <- vapply(result_out$Geneset.Peak.Genes, convert_ids, character(1))
      }
    }

    out_tsv <- args$output_tsv
    gs_tag <- sanitize_tag(gs)
    out_tsv <- sub("\\.tsv$", paste0(".", gs_tag, ".", m, ".tsv"), out_tsv)
    if (!grepl("\\.tsv$", out_tsv)) {
      out_tsv <- paste0(out_tsv, ".", gs_tag, ".", m, ".tsv")
    }

    write.table(result_out, file = out_tsv, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}
