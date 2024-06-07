# Add peak statistics to HOMER annotations file

library(readr)
library(dplyr)
library(tidyr)

args  <- commandArgs(trailingOnly=TRUE)

peaks_file  <- args[1]
homer_file  <- args[2]
combined_file  <- args[3]

# Load peak file
# the file extension is '.xls' but it's actually just a tsv file
peaks  <- read_tsv(peaks_file, comment = '#') %>%
    select(c("name","fold_enrichment","-log10(qvalue)"))

# Load HOMER annotations
homer  <- read_tsv(homer_file) %>%
    # rename first column to match peaks file
    rename(name = 1)

# Combine tables and write out
combined <- full_join(homer, peaks, by = 'name')
write_tsv(combined, file = combined_file)
