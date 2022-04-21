#!/usr/bin/env Rscript
rm(list=ls())

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ggvenn"))

# debug="TRUE"
# path="~/Documents/Projects/ccbr1155/CS030586/contrasts/siSmyd3_2m_Smyd3_0.25HCHO_500K_vs_siNC_2m_Smyd3_0.25HCHO_500K__dedup__narrowPeak"

# files=list.files(path=path)
debug="FALSE"
# path="/data/CCBR/projects/ccbr1155/CS030586_CARAP/diff"



if (debug){
  auc = grep(pattern = "_AUCbased_diffresults.txt",files)[1]
  aucresults = paste(path,files[auc],sep="/")
  frag = grep(pattern = "_fragmentsbased_diffresults.txt",files)[1]
  fragmentsresults = paste(path,files[frag],sep="/")
} else {
  # create parser object
  parser <- ArgumentParser()
  
  
  parser$add_argument("--aucresults", type="character", required=TRUE,
                      help="path to aucresults")
  parser$add_argument("--fragmentsresults", type="character", required=TRUE,
                      help="path to fragmentsresults")
  parser$add_argument("--pdf", type="character", required=TRUE,
                      help="output pdf file")
  parser$add_argument("--title", type="character", required=TRUE,
                      help="plot title in pdf")
  parser$add_argument("--fdr_cutoff", type="double", default=0.05, required=FALSE,
                      help="FDR cutoff [default %(default)s]")
  parser$add_argument("--log2fc_cutoff", type="double", default=0.59, required=FALSE,
                      help="log2foldchange cutoff [default %(default)s]")
  
  args <- parser$parse_args()
  aucresults=args$aucresults
  fragmentsresults=args$fragmentsresults
}

readandfilter <- function(resultsfile){
  df = read.csv(resultsfile,
                header = TRUE,
                row.names = NULL,
                sep = "\t",
                check.names = FALSE,
                comment.char = "#",
                strip.white = TRUE)
  df[is.na(df)] <- 1
  df$significance = "NS"
  k_fdr = df$padj < args$fdr_cutoff
  k_up = df$log2FoldChange > args$log2fc_cutoff
  k_down = df$log2FoldChange < (-1 * args$log2fc_cutoff)
  if (nrow(df[k_fdr & k_up,]) != 0) {df[k_fdr & k_up,]$significance = "UP"}
  if (nrow(df[k_fdr & k_down,]) != 0) {df[k_fdr & k_down,]$significance = "DOWN"}
  return(df)
}

auc_df = readandfilter(aucresults)
frag_df = readandfilter(fragmentsresults)

AUC_UP=auc_df[auc_df$significance=="UP",]$peakID
AUC_DOWN=auc_df[auc_df$significance=="DOWN",]$peakID
FRAG_UP=frag_df[frag_df$significance=="UP",]$peakID
FRAG_DOWN=frag_df[frag_df$significance=="DOWN",]$peakID
x<-list(AUC_UP=AUC_UP,AUC_DOWN=AUC_DOWN,FRAG_UP=FRAG_UP,FRAG_DOWN=FRAG_DOWN)
pdf(args$pdf)
ggvenn(x,
       stroke_size = 0.4, 
       set_name_size = 3,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")) +
  ggtitle(args$title) +
  theme(plot.title = element_text(size = 8))
dev.off()
