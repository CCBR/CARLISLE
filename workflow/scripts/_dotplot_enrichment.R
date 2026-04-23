#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

pick_col <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% names(df)) {
      return(nm)
    }
  }
  NULL
}

wrap_text <- function(x, width = 38) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
}

save_empty_plot <- function(outfile, title = "No GO Enrichment Results", subtitle = "Input TSV is empty.") {
  p_empty <- ggplot() +
    theme_void(base_size = 14) +
    annotate("text", x = 0.5, y = 0.58, label = title, size = 6, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.45, label = subtitle, size = 4.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  ggsave(outfile, plot = p_empty, width = 10, height = 7, dpi = 300, bg = "white")
  cat("Empty plot saved to", outfile, "\n")
}

save_no_data_and_quit <- function(outfile, subtitle) {
  save_empty_plot(outfile, subtitle = subtitle)
  quit(save = "no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript plot_enrichment_ggplot.R --input in.tsv --output out.png [--top_n 20] [--y_fontsize 8] [--wrap_width 38]")
}

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default)
  }
  if (idx[length(idx)] == length(args)) stop(paste0("Missing value for ", flag))
  args[idx[length(idx)] + 1]
}

infile <- get_arg("--input")
outfile <- get_arg("--output")
top_n <- as.integer(get_arg("--top_n", "20"))
y_fontsize <- as.numeric(get_arg("--y_fontsize", "8"))
wrap_width <- as.integer(get_arg("--wrap_width", "38"))

if (is.null(infile) || is.null(outfile)) {
  stop("Required: --input and --output")
}

if (!file.exists(infile)) {
  stop(paste("Input TSV does not exist:", infile))
}

infile_size <- file.info(infile)$size
if (is.na(infile_size) || infile_size == 0) {
  save_no_data_and_quit(outfile, "Input TSV is empty.")
}

df <- tryCatch(
  read.delim(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) {
    if (grepl("no lines available in input", e$message, fixed = TRUE)) {
      return(data.frame())
    }
    stop(paste("Failed to read input:", e$message))
  }
)

if (nrow(df) == 0) {
  save_no_data_and_quit(outfile, "No rows available in enrichment TSV.")
}

status_col <- pick_col(df, c("Status.x", "Status", "Status.Hybrid", "status"))
pvalue_col <- pick_col(df, c("P.value.x", "P.value", "P.value.Hybrid", "pvalue"))
fdr_col <- pick_col(df, c("FDR.Hybrid", "FDR", "adjPval", "padj", "qvalue"))
geneset_col <- pick_col(df, c("N.Geneset.Genes", "NumGenes", "GeneCount", "Count"))
peak_col <- pick_col(df, c("N.Geneset.Peak.Genes", "NumPeaks", "num.peaks", "Count"))
ratio_col <- pick_col(df, c("GeneRatio", "gene_ratio"))
desc_col <- pick_col(df, c("Description", "Term", "Pathway", "Name", "Geneset.ID"))

missing <- c()
if (is.null(status_col)) missing <- c(missing, "status column")
if (is.null(pvalue_col)) missing <- c(missing, "p-value column")
if (is.null(geneset_col)) missing <- c(missing, "geneset size column")
if (is.null(peak_col) && is.null(ratio_col)) missing <- c(missing, "peak count or GeneRatio column")
if (is.null(desc_col)) missing <- c(missing, "description column")
if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

df <- df[tolower(trimws(as.character(df[[status_col]]))) == "enriched", , drop = FALSE]
if (nrow(df) == 0) save_no_data_and_quit(outfile, "No enriched pathways found.")

df[[pvalue_col]] <- suppressWarnings(as.numeric(df[[pvalue_col]]))
df <- df[is.finite(df[[pvalue_col]]) & df[[pvalue_col]] > 0, , drop = FALSE]
if (nrow(df) == 0) save_no_data_and_quit(outfile, "No enriched pathways with valid p-values.")

if (is.null(fdr_col)) fdr_col <- pvalue_col

df[[fdr_col]] <- suppressWarnings(as.numeric(df[[fdr_col]]))
df[[geneset_col]] <- suppressWarnings(as.numeric(df[[geneset_col]]))
if (!is.null(peak_col)) df[[peak_col]] <- suppressWarnings(as.numeric(df[[peak_col]]))
if (!is.null(ratio_col)) df[[ratio_col]] <- suppressWarnings(as.numeric(df[[ratio_col]]))

df <- df[is.finite(df[[geneset_col]]) & df[[geneset_col]] > 0, , drop = FALSE]
if (nrow(df) == 0) save_no_data_and_quit(outfile, "No rows with positive geneset size.")

if (!is.null(ratio_col)) {
  df$GeneRatio <- df[[ratio_col]]
} else {
  df$GeneRatio <- df[[peak_col]] / df[[geneset_col]]
}

if (!is.null(peak_col)) {
  df$GeneCount <- df[[peak_col]]
} else {
  df$GeneCount <- df[[geneset_col]]
}

df <- df[
  is.finite(df$GeneRatio) & df$GeneRatio > 0 &
    is.finite(df$GeneCount) & df$GeneCount > 0 &
    is.finite(df[[fdr_col]]), ,
  drop = FALSE
]
if (nrow(df) == 0) save_no_data_and_quit(outfile, "No plottable rows after filtering.")

fdr_vals <- suppressWarnings(as.numeric(df[[fdr_col]]))
p_vals <- suppressWarnings(as.numeric(df[[pvalue_col]]))

# Keep log10-safe FDR values; fall back to p-values when FDR is non-positive.
fdr_plot <- ifelse(is.finite(fdr_vals) & fdr_vals > 0, fdr_vals, NA_real_)
fdr_plot <- ifelse(is.na(fdr_plot) & is.finite(p_vals) & p_vals > 0, p_vals, fdr_plot)
if (all(!is.finite(fdr_plot))) {
  save_no_data_and_quit(outfile, "No positive FDR/p-value values available for plotting.")
}

min_positive <- min(fdr_plot[is.finite(fdr_plot)], na.rm = TRUE)
fdr_plot[!is.finite(fdr_plot)] <- min_positive
df$FDR_plot <- fdr_plot

df <- df[order(df[[fdr_col]], decreasing = FALSE), , drop = FALSE]
df <- head(df, top_n)
df <- df[order(df$GeneRatio, decreasing = TRUE), , drop = FALSE]

df$Term <- as.character(df[[desc_col]])
df$TermWrapped <- wrap_text(df$Term, width = wrap_width)
# Wrapped term labels can collide across rows; deduplicate factor levels while
# preserving display order so repeated pathway names do not abort plotting.
df$TermWrapped <- factor(df$TermWrapped, levels = rev(unique(df$TermWrapped)))
df$FDR <- df$FDR_plot

fdr_min <- min(df$FDR, na.rm = TRUE)
fdr_max <- max(df$FDR, na.rm = TRUE)
exp_min <- floor(log10(fdr_min))
exp_max <- ceiling(log10(fdr_max))
legend_breaks <- 10^(seq(exp_min, exp_max))
legend_labels <- paste0("1e", formatC(log10(legend_breaks), format = "d"))

fig_h <- max(7, 0.35 * nrow(df))

p <- ggplot(df, aes(x = GeneRatio, y = TermWrapped)) +
  geom_point(aes(size = GeneCount, color = FDR), shape = 16, alpha = 0.95) +
  scale_color_gradientn(
    colours = c("#b2182b", "#f7f7f7", "#2166ac"),
    trans = "log10",
    name = "FDR",
    limits = c(10^exp_min, 10^exp_max),
    breaks = legend_breaks,
    labels = legend_labels
  ) +
  scale_size_continuous(name = "Gene Count", range = c(2.5, 12)) +
  labs(x = "GeneRatio", y = NULL) +
  guides(
    color = guide_colorbar(order = 1, reverse = TRUE),
    size = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.4),
    panel.grid.minor = element_line(color = "#F7F7F7", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = y_fontsize),
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = c(0, 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(outfile, plot = p, width = 10, height = fig_h, dpi = 300, bg = "white")
cat("Plot saved to", outfile, "\n")
