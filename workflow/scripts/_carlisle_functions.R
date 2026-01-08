########################################################################
# LIBRARY
########################################################################
library(tidyverse)
load_packages <- function() {
  pkgs <- "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
chipenrich
ChIPseeker
DESeq2
edgeR
ELBOW
EnhancedVolcano
GenomicFeatures
ggfortify
ggrepel
htmltools
HTSFilter
latticeExtra
org.Hs.eg.db
org.Mm.eg.db
pander
pdp
plotly
RColorBrewer
reshape2
rtracklayer
tidyverse
TxDb.Hsapiens.UCSC.hg19.knownGene
TxDb.Hsapiens.UCSC.hg38.knownGene
TxDb.Mmusculus.UCSC.mm10.knownGene
yaml"
  package_vctr <- pkgs %>%
    stringr::str_split("\n") %>%
    unlist()
  invisible(lapply(package_vctr, library, character.only = TRUE))
}

########################################################################
# SPIKE IN PLOT
########################################################################
GENERATE_SPIKEIN_PLOT <- function(input_df, spike_type) {
  spike_df <- data.frame()
  for (rowid in rownames(input_df, spike_type)) {
    # spike_type="NC_000913.3"

    # read in file
    stats <- read.table(input_df[rowid, "bam"])
    stats <- stats[, c("V1", "V3")]
    colnames(stats) <- c("location", "read_count")

    # add metadata
    stats$sampleid <- input_df[rowid, "repid"]
    stats$groupid <- input_df[rowid, "sampleid"]

    if (nrow(stats) == 0) {
      spike_df <- subset(stats, location == spike_type)
    } else {
      spike_df <- rbind(
        subset(stats, location == spike_type),
        spike_df
      )
    }
  }

  p <- ggplot(data = spike_df, aes(x = sampleid, y = read_count, fill = groupid)) +
    geom_bar(stat = "identity")
  p_final <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste0("Spike-in control values\n", spike_type))
  print(p_final)

  return(spike_df)
}

########################################################################
# GO ENRICHMENT
########################################################################
READ_PEAK_FILE <- function(peak_file_in) {
  # `peak_file_in` may be a 1x1 tibble/data.frame cell (e.g. peak_df[rowid, "peak_bed"]).
  # Normalize to a single filepath string before reading.
  if (is.data.frame(peak_file_in)) {
    peak_file_in <- peak_file_in[[1]]
  }
  peak_file_in <- as.character(peak_file_in)[1]

  peak_df <- read.csv(peak_file_in, sep = "\t", header = FALSE)[, c("V1", "V2", "V3")]
  colnames(peak_df) <- c("chrom", "start", "end")

  return(peak_df)
}

SET_LOCUST_DEF <- function(locus_loc_in) {
  # depending on where the majority of peaks fall, ID locus choice to use
  if (locus_loc_in == "<0.1" | locus_loc_in == "0.1 - 1") {
    locus_loc_in <- "1kb"
  } else if (locus_loc_in == "1 - 5") {
    locus_loc_in <- "5kb"
  } else if (locus_loc_in == "5 - 10") {
    locus_loc_in <- "10kb"
  } else {
    locus_loc_in <- "not applicable given distribution"
  }

  print(paste0("The locust defintion is determined to be: ", locus_loc_in))
  return(locus_loc_in)
}

SET_LOCUST_LIST <- function(locus_loc_in) {
  # depending on where the majority of peaks fall, ID locus choice to use
  if (locus_loc_in == "<0.1" | locus_loc_in == "0.1 - 1") {
    locusdf_list <- c("1kb", "1kb_outside", "1kb_outside_upstream")
  } else if (locus_loc_in == "1 - 5") {
    locusdf_list <- c("5kb", "5kb_outside", "5kb_outside_upstream")
  } else if (locus_loc_in == "5 - 10") {
    locusdf_list <- c("10kb", "10kb_outside", "10kb_outside_upstream")
  } else {
    locusdf_list <- "none"
  }

  return(locusdf_list)
}

PLOT_QC_FUNCTIONS <- function(function_in, rowid_in, l_id) {
  if (function_in == "polyenrich") {
    p <- plot_polyenrich_spline(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid_in]),
      locusdef = l_id, genome = speciesID
    )
  } else if (function_in == "spline") {
    p <- plot_chipenrich_spline(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid_in]),
      locusdef = l_id, genome = speciesID
    )
  } else if (function_in == "cov") {
    p <- plot_gene_coverage(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid_in]),
      locusdef = l_id, genome = speciesID
    )
  } else {
    print(paste0("Missing function", function_in))
  }
  return(p)
}

PLOT_QC_MAIN <- function(function_in, rowid_in) {
  locusdf_list_str <- as.character(peak_df[["locusdf_list"]][1])[1]
  locusdf_list <- c(base::strsplit(locusdf_list_str, ",")[[1]], "nearest_tss", "nearest_gene")
  locusdf_list <- locusdf_list[locusdf_list != "none"]

  # create plots, catching for errors
  counter <- 1
  legend_text <- ""
  row_count <- 1
  col_count <- 1
  p <- list()
  for (l_id in locusdf_list) {
    plot_test <- tryCatch(
      {
        PLOT_QC_FUNCTIONS(function_in, rowid_in, l_id)
      },
      warning = function(w) {
        "warning"
      },
      error = function(e) {
        "failed"
      }
    )

    # save successful plots in plot_list
    # save successful plots legend
    if (class(plot_test) == "trellis") {
      p[[counter]] <- plot_test
      counter <- counter + 1

      # if true, this is a new row
      if (col_count == 1) {
        legend_text <- paste0(legend_text, "Row ", row_count, ":  | Col 1 (", l_id, ")")
        col_count <- 2
      } else if (col_count == 2 & row_count == 1) {
        legend_text <- paste0(legend_text, " | Col 2 (", l_id, ")\n")
        row_count <- row_count + 1
        col_count <- 1
      } else if (col_count == 2 & row_count > 1) {
        legend_text <- paste0(legend_text, " | Col 2 (", l_id, ")")
        col_count <- col_count + 1
      } else {
        legend_text <- paste0(legend_text, " | Col 3 (", l_id, ")\n")
        col_count <- 1
        row_count <- row_count + 1
      }
    }
  }

  # print with data key
  cat(legend_text)
  n <- length(p)
  if (n == 0) {
    print("All plots failed for this sample")
  } else if (n == 1) {
    plot(c(p[[1]]))
  } else if (n == 2) {
    plot(c(p[[1]], p[[2]]))
  } else if (n == 3) {
    plot(c(p[[1]], p[[2]], p[[3]]))
  } else if (n == 4) {
    plot(c(p[[1]], p[[2]], p[[3]], p[[4]]))
  } else if (n == 5) {
    plot(c(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]]))
  } else if (n == 6) {
    plot(c(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]]))
  } else {
    plot(paste0("Missing N", n))
  }
}

GO_ANALYSIS_MAIN <- function(rowid, peak_enrichment) {
  sampleid <- peak_df[["sampleid"]][rowid]
  peakcaller <- peak_df[["peak_caller"]][rowid]
  peaktype <- peak_df[["peak_type"]][rowid]
  locus_def <- peak_df[["locus_loc_short"]][rowid]

  n_cores_in <- 1
  if (exists("carlisle_n_cores", inherits = TRUE)) {
    n_cores_in <- as.integer(get("carlisle_n_cores", inherits = TRUE))
  }

  print(paste0("** ", sampleid, " | ", peakcaller, " | ", peaktype, " **"))
  print(paste0("-- Peak analysis: ", peak_enrichment))

  if (locus_def == "none") {
    locus_def <- "nearest_tss"
  }

  if (peak_enrichment == "hybrid") {
    results <- hybridenrich(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      genome = speciesID, genesets = geneset_id,
      locusdef = locus_def, qc_plots = F, # randomization = 'complete',
      out_name = NULL, n_cores = n_cores_in, min_geneset_size = 10
    )
  } else if (peak_enrichment == "broad") {
    results <- broadenrich(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      genome = speciesID, genesets = geneset_id,
      locusdef = locus_def, qc_plots = FALSE, # randomization = 'complete',
      out_name = NULL, n_cores = n_cores_in, min_geneset_size = 10
    )
  } else if (peak_enrichment == "enrich") {
    results <- chipenrich(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      genome = speciesID, genesets = geneset_id,
      locusdef = locus_def, qc_plots = FALSE, # randomization = 'complete',
      out_name = NULL, n_cores = n_cores_in
    )
  } else if (peak_enrichment == "poly") {
    results <- polyenrich(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      genome = speciesID, genesets = geneset_id,
      method = "polyenrich", locusdef = locus_def,
      qc_plots = FALSE, out_name = NULL, n_cores = n_cores_in
    )
  } else if (peak_enrichment == "poly_weighted") {
    results <- polyenrich(
      peaks = READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      genome = speciesID, genesets = geneset_id,
      method = "polyenrich_weighted", locusdef = "enhancer_plus5kb",
      qc_plots = FALSE, out_name = NULL, n_cores = n_cores_in
    )
  } else if (peak_enrichment == "reglocation") {
    results <- proxReg(READ_PEAK_FILE(peak_df[["peak_bed"]][rowid]),
      reglocation = "tss", genome = speciesID,
      genesets = geneset_id, out_name = NULL
    )
  }

  # set results, sig
  result.out <- results$results
  alpha <- sum(result.out$P.value < 0.05) / nrow(result.out)
  print(paste0("--sig: ", alpha))

  # write out (TSV for consistency)
  fpath <- paste0(output_dir, "/", sampleid, ".", peakcaller, ".", dedup_status, ".", peaktype, ".", peak_enrichment, "_", geneset_id, ".tsv")
  write.table(result.out, file = fpath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
