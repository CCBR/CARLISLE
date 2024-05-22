# Add peak statistics to HOMER annotations file

library(openxlsx)

args = commandArgs(trailingOnly=TRUE)

peaks.file = args[1]
homer.file = args[2]
combined.file = args[3]

# Load peak file
peaks = read.table(peaks.file,skip = 23,header=TRUE,sep='\t',check.names = FALSE)

# Load HOMER annotations
homer = read.table(homer.file,header=TRUE,sep='\t',check.names=FALSE,comment.char="",quote="")

# Combine tables and write out
combined = merge(homer,peaks[,c("name","fold_enrichment","-log10(qvalue)")],by.x=1,by.y="name")
write.xlsx(combined,file=combined.file)
