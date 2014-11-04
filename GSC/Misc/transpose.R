setwd("C:/Documents and Settings/obig/Desktop/Trevor")
#datafile="coverage_adjusted_mincov2_gene_summary.txt"
datafile="preRx_coding_mutations_090615IB.txt"

rawdata=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
#data_subset=rawdata[,c(1,9:length(colnames(rawdata)))]
data_subset=rawdata


data=t(data_subset)
#write.table(data, file="coverage_adjusted_mincov2_gene_summary.RFclean.txt", row.names=TRUE, col.names=FALSE,quote=FALSE, sep="\t")
write.table(data, file="preRx_coding_mutations_090615IB.RFclean.txt", row.names=TRUE, col.names=FALSE,quote=FALSE, sep="\t")

