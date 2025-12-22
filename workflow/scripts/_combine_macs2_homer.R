# Add peak statistics to HOMER annotations file

library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

peaks_file <- args[1]
homer_file <- args[2]
combined_tsv <- args[3]
combined_xlsx <- args[4]

# Load peak file
# the file extension is '.xls' but it's actually just a tsv file
peaks <- read_tsv(peaks_file, comment = "#")

# Check if peaks file has data
if (nrow(peaks) > 0 && "name" %in% colnames(peaks)) {
  peaks <- peaks %>%
    select(c("name", "fold_enrichment", "-log10(qvalue)"))
} else {
  # Create empty dataframe with expected columns
  peaks <- data.frame(
    name = character(0),
    fold_enrichment = numeric(0),
    `-log10(qvalue)` = numeric(0),
    check.names = FALSE
  )
}

# Load HOMER annotations
homer <- read_tsv(homer_file)

# Check if HOMER file has data
if (nrow(homer) > 0 && ncol(homer) > 0) {
  # rename first column to match peaks file
  homer <- homer %>% rename(name = 1)
} else {
  # Create empty dataframe with just the name column
  homer <- data.frame(name = character(0))
}

# Combine tables and write out
combined <- full_join(homer, peaks, by = "name")
write_tsv(combined, file = combined_tsv)
openxlsx::write.xlsx(combined, file = combined_xlsx)
