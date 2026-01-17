# Plot HOMER enrichment

library(reshape2)
library(ggplot2)
library(plyr)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

peaks.dir <- args[1]
peak_mode <- args[2]
dupstatus <- args[3]
fig <- args[4]

# Load HOMER annotation summary
homer.files <- list.files(
  path = peaks.dir, pattern = "annotation.summary",
  recursive = TRUE, full.names = TRUE
)

# Read annotation files, handling empty files
homer <- lapply(homer.files, function(x) {
  df <- read.table(x, sep = "\t", header = TRUE, check.names = FALSE, nrows = 14)
  # Return NULL if no data rows (only header)
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
})

# Remove NULL entries (files with no peaks)
names(homer) <- gsub(
  "gopeaks/|macs2/|seacr/", "",
  gsub(
    ".annotation.summary", "",
    gsub(
      "/annotation/homer", "",
      gsub("^.*results/peaks/", "", homer.files)
    )
  )
)
homer <- homer[!sapply(homer, is.null)]

# Check if we have any valid data
if (length(homer) == 0) {
  warning("No valid annotation data found. Creating empty plot.")
  png(fig, width = 700, height = 400)
  plot.new()
  text(0.5, 0.5, "No peaks with HOMER annotation available", cex = 1.5)
  dev.off()
  quit(save = "no", status = 0)
}

homer <- lapply(homer, function(x) x[1:rownames(x[which(x$Annotation == "Annotation"), ]) - 1, ])
homer <- ldply(homer, .id = "File")

homer <- cbind(homer, colsplit(homer$File, "/", c("Threshold", "Sample")))
homer <- cbind(homer, colsplit(homer$Sample, "[.]", c("Comparison", "Duplication", "Caller")))
homer <- cbind(homer, colsplit(homer$Comparison, "_vs_", c("Replicate", "Control")))
homer <- homer[which(homer$Duplication == dupstatus & homer$Caller == peak_mode), ]
homer$`LogP enrichment (+values depleted)` <- as.numeric(homer$`LogP enrichment (+values depleted)`)

# Plot enrichment heatmap
homer <- dcast(homer, Replicate ~ Annotation, value.var = "LogP enrichment (+values depleted)")
rownames(homer) <- homer$Replicate
homer <- homer[, 2:dim(homer)[2]]

png(fig, width = 700, height = 400)
hm <- Heatmap(homer,
  name = "LogP enrichment\n(+values depleted)",
  col = colorRamp2(
    breaks = c(-70000, -50, 0, 50, 70000),
    colors = c("red", "yellow", "white", "green", "blue")
  ),
  width = unit(100, "mm")
)
draw(hm, heatmap_legend_side = "left")
dev.off()
