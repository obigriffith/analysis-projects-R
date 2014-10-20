setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/SUM149")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/A431")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/BT474")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/HS578T")


datafile="SUM149_compare_mode_3_standard.peaks"
datafile="A431_compare_mode_3_standard.peaks"
datafile="BT474_compare_mode_3_standard.peaks"
datafile="HS578T_compare_mode_3_standard.peaks"

#data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", row.names=NULL)

pdf("SUM149_peak_height_dist.pdf")
pdf("A431_peak_height_dist.pdf")
pdf("BT474_peak_height_dist.pdf")
pdf("HS578T_peak_height_dist.pdf")

#summary=table(data[,5])
summary=table(data[,6])
barplot(summary, cex.names=0.45, col="blue", xlab="peak height", ylab="frequency", main="SUM149 iressa-treated me1 peak height distribution")
barplot(summary, cex.names=0.45, col="blue", xlab="peak height", ylab="frequency", main="A431 iressa-treated me1 peak height distribution")
barplot(summary, cex.names=0.45, col="blue", xlab="peak height", ylab="frequency", main="BT474 iressa-treated me1 peak height distribution")
barplot(summary, cex.names=0.45, col="blue", xlab="peak height", ylab="frequency", main="HS578T iressa-treated me1 peak height distribution")
dev.off()


pdf("SUM149_peak_height_dist_log10.pdf")
pdf("A431_peak_height_dist_log10.pdf")
pdf("BT474_peak_height_dist_log10.pdf")
pdf("HS578T_peak_height_dist_log10.pdf")

barplot(log10(summary), cex.names=0.45, col="blue", xlab="peak height", ylab="log10 frequency", main="SUM149 iressa-treated me1 peak height distribution")
barplot(log10(summary), cex.names=0.45, col="blue", xlab="peak height", ylab="log10 frequency", main="A431 iressa-treated me1 peak height distribution")
barplot(log10(summary), cex.names=0.45, col="blue", xlab="peak height", ylab="log10 frequency", main="BT474 iressa-treated me1 peak height distribution")
barplot(log10(summary), cex.names=0.45, col="blue", xlab="peak height", ylab="log10 frequency", main="HS578T iressa-treated me1 peak height distribution")
dev.off()
