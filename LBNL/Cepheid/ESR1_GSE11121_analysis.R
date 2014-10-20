library("gplots")
library("heatmap.plus")
library(genefilter)
library(gplots)

#datafile="C:/Users/Obi/Documents/Projects/Cepheid/processed/customCDF/GSE11121_gcrma.txt"
datafile="C:/Users/Obi/Documents/Projects/Cepheid/processed/standardCDF/GSE11121_gcrma.txt"
clindatafile="C:/Users/Obi/Documents/Projects/Cepheid/filtered_GSE11121_data_anno.txt"

#outdir="C:/Users/Obi/Documents/Projects/Cepheid/analysis/customCDF/GSE11121_samples"
outdir="C:/Users/Obi/Documents/Projects/Cepheid/analysis/standardCDF/GSE11121_samples"
outfile3B = "GSE11121_gcrma_heatmap_ESR1_ordered.pdf"
outfile4B = "GSE11121_gcrma_waterfall_ESR1.pdf"
outfile5B = "GSE11121_gcrma_hist_ESR1.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clin_data=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
raw_data=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

#Change to output dir
setwd(outdir)
header=colnames(raw_data)

#Preliminary gene filtering
X=raw_data[,4:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=raw_data[filt,] 

#Single-gene Data
ESR1_Data=raw_data[which(raw_data[,3]=="ESR1"),]

#Heatmap (single color sidebar) - ESR1 only, order samples by decreasing ESR1 instead of clustering. Color side bar requires same reorder
x=sort(as.matrix(ESR1_Data[,4:length(header)]), decreasing=TRUE)
z=order(as.matrix(ESR1_Data[,4:length(header)]), decreasing=TRUE)
x=rbind(x,x) #duplicate to trick heatmap.2 which requires minimum of 2 x 2 matrix
pdf(file=outfile3B)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE, labRow=FALSE, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Waterfall plot - ESR1
y=sort(as.numeric(ESR1_Data[,4:length(header)]), decreasing=TRUE)
x=0:(length(y)-1)
pdf(file=outfile4B)
plot(y, type="n", main="ESR1 waterfall for all GSE11121 samples", ylab="Log2 GCRMA value", xlab="Sample")
polygon(c(min(x), x, max(x), 0), c(0, y, 0, 0), col="blue")
dev.off()

#Histogram - ESR1
pdf(file=outfile5B)
hist(as.numeric(ESR1_Data[,4:length(header)]), breaks=40, col="blue", main="ESR1 histogram for GSE11121 samples", xlab="Log2 GCRMA value", ylab="Frequency")
dev.off()

