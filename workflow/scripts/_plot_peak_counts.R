# Plot number of peaks called per sample by caller, duplication status, and threshold

library(reshape2)
library(ggplot2)
library(plyr)
library(openxlsx)

theme_set(theme_bw())

args = commandArgs(trailingOnly=TRUE)

input = args[1]

# Load number of peaks
print("Load peaks")
peaks = read.table(input,sep="")
colnames(peaks) = c("Peaks","File")
peaks = peaks[which(peaks$File != "total"),]
peaks$File = gsub(".peaks.bed$","",gsub("/peak_output","",gsub("^.*results/peaks/","",peaks$File)))
peaks = cbind(peaks,colsplit(peaks$File,"/",c("Threshold","Caller","Sample")))
peaks = cbind(peaks,colsplit(peaks$Sample,"[.]",c("Comparison","Duplication","Mode")))
peaks$Caller = paste(peaks$Caller,peaks$Mode,sep="_")
peaks = peaks[,c("Comparison","Caller","Duplication","Threshold","Peaks")]
peaks = cbind(peaks,colsplit(peaks$Comparison,"_vs_",c("Replicate","Control")))
peaks$Sample = gsub("_[1-9]$","",peaks$Replicate)

setwd(args[2])

# Write out table
print("Write table")
peak_output = dcast(peaks,Threshold+Duplication+Comparison~Caller,value.var="Peaks")
peak_output = peak_output[order(peak_output$Threshold,peak_output$Duplication,peak_output$Comparison),]
write.xlsx(peak_output,file="Peak counts.xlsx")

# For each threshold, plot number of peaks
for (threshold in unique(peaks$Threshold)){
	print(threshold)
peaks_threshold = peaks[which(peaks$Threshold == threshold),]

# Compare peaks with and without duplication
peaks_threshold_dup = dcast(peaks_threshold,Replicate+Caller~Duplication,value.var="Peaks")
cor=round(as.numeric(unlist(cor.test(peaks_threshold_dup$dedup,peaks_threshold_dup$no_dedup))["estimate.cor"]),2)
peaks.max=max(max(peaks_threshold_dup$dedup),max(peaks_threshold_dup$no_dedup))

png(paste("duplication_corr.",threshold,".png",sep=""),width = 350, height = 300)
    print(ggplot(peaks_threshold_dup,aes(x=log10(no_dedup),y=log10(dedup))) + geom_point(aes(color=Caller)) +
        geom_abline(intercept=0,slope=1,linetype="dashed") +
        xlab("log10(Peaks), no deduplication") + ylab("log10(Peaks), deduplication") +
        xlim(0,log10(peaks.max)+0.5) + ylim(0,log10(peaks.max)+0.5) +
        ggtitle(paste("Peaks by read duplication, q-value threshold ",threshold,sep=""),
	    subtitle=paste("Pearson correlation: ",cor,sep="")))
dev.off()

# Plot summary across callers
width=(length(unique(peaks_threshold$Sample))*50)+100
height=(length(unique(peaks_threshold$Caller))*100)+100
png(paste("peaks_by_caller.",threshold,".png",sep=""),width = width, height = height)
    print(ggplot(peaks_threshold[which(peaks_threshold$Duplication == "dedup"),],aes(x=Sample,y=Peaks)) +
  	geom_boxplot() +
	facet_wrap(~Caller,nrow=length(unique(peaks_threshold$Caller)),scales="free_y") +
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
 	ggtitle(paste("Peaks,\nq-value threshold ",threshold,",\ndeduplicated reads",sep="")))	
dev.off()
}
