# Plot correlation heatmap and PCA

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

args = commandArgs(trailingOnly=TRUE)

in.corr = args[1]
in.pca = args[2]
in.read_depth = args[3]
dupstatus = args[4]
out.corr = args[5]
out.pca = args[6]

# Read deeptools output
print("Load heatmap")
heatmap = read.table(in.corr,check.names = FALSE)
colnames(heatmap) = gsub("[.]no_dedup","",gsub("[.]dedup","",colnames(heatmap)))
rownames(heatmap) = gsub("[.]no_dedup","",gsub("[.]dedup","",rownames(heatmap)))
print(heatmap)

print("Load PCA")
pca = read.table(in.pca,header=TRUE,check.names = FALSE)
print(pca)

# Create metadata table
print("Create metadata table")
metadata = data.frame(Replicate=colnames(heatmap))
metadata$Sample = gsub("_[1-9]$","",metadata$Replicate)

# Add read depth
read_depth = read.table(in.read_depth,sep='\t',header=TRUE)
metadata = merge(metadata,read_depth[,c("sample_name","nreads","no_dedup_nreads_genome","dedup_nreads_genome")],
                 by.x="Replicate",by.y="sample_name")
if(dupstatus == "dedup"){
      metadata$Depth = metadata$dedup_nreads_genome/1000000
}else{
      metadata$Depth = metadata$no_dedup_nreads_genome/1000000
}
metadata = metadata[match(rownames(heatmap),metadata$Replicate),]
print(metadata)

# Create colors
sample_colors = setNames(colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(metadata$Sample))),
                       unique(metadata$Sample))
depth_colors = colorRamp2(breaks=c(0,30,60),colors=c("white","grey50","black"))

colors = list(Sample=sample_colors,
              Depth=depth_colors)

# Heatmap

## Plot heatmap
print("Print heatmap")
png(out.corr,width=600,height=600)
print(Heatmap(heatmap,
        top_annotation = HeatmapAnnotation(df=metadata[,c("Sample","Depth")],
                                           col=colors[c("Sample","Depth")]),
        show_row_dend = FALSE,
        show_row_names = FALSE,
        col = colorRamp2(c(0,1), c("white", "blue")),
	heatmap_legend_param = list(
          title = "Pearson\ncorrelation"
        ),
        column_dend_height = unit(0.8,"inches"),
        name="Pearson correlation, genome-wide coverage (10kb bins)"))
dev.off()

# PCA
print("Print PCA")

## Calculate variance from eigenvalue
eigenvalue = pca[,c("Component","Eigenvalue")]
eigenvalue$Variance = eigenvalue$Eigenvalue/sum(eigenvalue$Eigenvalue)*100

pca = melt(pca[,1:dim(pca)[2]-1],id.var="Component",variable.name="Replicate",value.name="Loading")
pca = dcast(pca,Replicate~Component,value.var="Loading")
pca$Replicate = gsub("[.]no_dedup","",gsub("[.]dedup","",pca$Replicate))
pca = merge(pca,metadata,by="Replicate")

## Plot PCA
png(out.pca)
print(ggplot(pca,aes(x=`1`,y=`2`,color=Sample)) + geom_point(size=3) +
  xlab(paste("PC1 (",round(eigenvalue$Variance[1],1),"%)",sep="")) +
  ylab(paste("PC2 (",round(eigenvalue$Variance[2],1),"%)",sep="")) +
  scale_color_manual(values=sample_colors) +
  ggtitle("PCA, genome-wide coverage (10kb bins)") +
  theme_classic() + theme(legend.position="bottom"))
dev.off()
