library("gplots")
library("heatmap.plus")
library(randomForest)
library(ROCR)

#Cell line info
celllinefile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/cell_line_info.2.txt"

#Gene
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"

#Transcript
#datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Matrix_TranscriptExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Expressed_TranscriptExpression_v53.txt"

#KnownJunction
#datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Matrix_KnownJunctionExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Expressed_KnownJunctionExpression_v53.txt"

#Set output dir for feature type
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/analysis/TFAP2C/"
setwd(outdir)

TFAP2Camp_heatmap_outfile = "BCCL_TFAP2CAmp_67libs_heatmap.pdf"
TFAP2C_heatmap_outfile = "BCCL_TFAP2C_67libs_heatmap.pdf"
waterfall_file = "BCCL_TFAP2C_67libs_waterfall.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_cell_data=read.table(celllinefile, header = TRUE, na.strings = "NA", sep="\t")
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:length(colnames(raw_data_import))]
raw_exp_status=raw_exp_status_import[,5:length(colnames(raw_data_import))]

#Make sure that cell line info and raw data are ordered the same!!!
libs=colnames(raw_data)
lib_names=as.vector(raw_cell_data[,"Sample.Name"])
cbind(lib_names,libs)

#If ok, use clean names for data object
colnames(raw_data)=lib_names

#All libraries
#data=raw_data
#feat_data=raw_feat_data
#exp_status=raw_exp_status
#cell_data=raw_cell_data

#Exclude low quality 
high_qual=which(as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int")

#Exclude low quality and Normal
#high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"subtype"]!="Normal"))

#Exclude low quality and "Unknown/Normal subtype" libs
#high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"subtype"]!="Unknown") & as.vector(raw_cell_data[,"subtype"]!="Normal"))

#Filter down to only libraries that also have exon Array data for comparisons purposes
#RNA_exon_libs=c("184A1","184B5","600MPE","BT20","BT474","BT483","BT549","CAMA1","DU4475","HCC1143","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1599","HCC1806","HCC1937","HCC1954","HCC202","HCC2218","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM225CWN","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")
#high_qual=which(lib_names %in% RNA_exon_libs)

#Apply library filter to datasets
data=raw_data[,high_qual]
feat_data=raw_feat_data
exp_status=raw_exp_status[,high_qual]
cell_data=raw_cell_data[high_qual,]
libs=colnames(data)

#Retrieve cell line details
lib_names=as.vector(cell_data[,"Sample.Name"])
subtypes=as.vector(cell_data[,"subtype"])
ERBB2=as.vector(cell_data[,"ERBB2New"])
Qualities=as.vector(cell_data[,"Quality"])
cbind(lib_names,libs)


#Set up colors for sidebars
colors=rainbow(6)
#Subtype
subtype_colors=subtypes
subtype_colors[subtype_colors=="Basal"]=colors[1]
subtype_colors[subtype_colors=="Luminal"]=colors[2]
subtype_colors[subtype_colors=="ClaudinLow"]=colors[3]
subtype_colors[subtype_colors=="Basal_NM"]=colors[4]
subtype_colors[subtype_colors=="Normal"]=colors[5]
subtype_colors[subtype_colors=="Unknown"]=colors[6]
names(subtype_colors)=lib_names #Set names so that it can be used with reordered cell line lists later

#TFAP2C status
her2_colors=ERBB2
her2_colors[her2_colors=="Amp"]=colors[1]
her2_colors[her2_colors=="NoAmp"]=colors[5]

#Quality status
quality_colors=Qualities
quality_colors[quality_colors=="High"]=colors[1]
quality_colors[quality_colors=="Int"]=colors[3]
quality_colors[quality_colors=="Low"]=colors[5]

#ERBB2 data
ERBB2_ind=which(feat_data[,"Seq_Name"]=="ERBB2")

#Get TFAP2C Amp genes (chose genes from UCSC browser, hg19, zoomed out, centered on TFAP2C)
CBLN4_ind=which(feat_data[,"Seq_Name"]=="CBLN4")
MC3R_ind=which(feat_data[,"Seq_Name"]=="MC3R")
C20orf108_ind=which(feat_data[,"Seq_Name"]=="C20orf108")
AURKA_ind=which(feat_data[,"Seq_Name"]=="AURKA")
CSTF1_ind=which(feat_data[,"Seq_Name"]=="CSTF1")
CASS4_ind=which(feat_data[,"Seq_Name"]=="CASS4")
C20orf43_ind=which(feat_data[,"Seq_Name"]=="C20orf43")
GCNT7_ind=which(feat_data[,"Seq_Name"]=="GCNT7")
C20orf106_ind=which(feat_data[,"Seq_Name"]=="C20orf106")
C20orf107_ind=which(feat_data[,"Seq_Name"]=="C20orf107")
TFAP2C_ind=which(feat_data[,"Seq_Name"]=="TFAP2C")
BMP7_ind=which(feat_data[,"Seq_Name"]=="BMP7")
SPO11_ind=which(feat_data[,"Seq_Name"]=="SPO11")
RAE1_ind=which(feat_data[,"Seq_Name"]=="RAE1")

AmpGenes_ind=c(CBLN4_ind,MC3R_ind,C20orf108_ind,AURKA_ind,CSTF1_ind,CASS4_ind,C20orf43_ind,GCNT7_ind,C20orf106_ind,C20orf107_ind,TFAP2C_ind,BMP7_ind,SPO11_ind,RAE1_ind)
TFAP2C_ERBB2_ind=c(TFAP2C_ind,ERBB2_ind)

data_TFAP2Camp=data[AmpGenes_ind,]
feat_data_TFAP2Camp=feat_data[AmpGenes_ind,]
exp_status_TFAP2Camp=exp_status[AmpGenes_ind,]
rownames(data_TFAP2Camp)=feat_data_TFAP2Camp[,"Seq_Name"]
gene_names=rownames(data_TFAP2Camp)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap - TFAP2C Amp
#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(data_TFAP2Camp+1)
pdf(file=TFAP2Camp_heatmap_outfile)
x=as.matrix(z)
main_title="TFAP2C Amplicon RNAseq gene-level expression"
clab=cbind(subtype_colors,her2_colors,quality_colors)
colnames(clab) = c("Subtype","Her2","Quality")
#heatmap.plus(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Rowv=NA, main=main_title, ColSideColors=clab, labRow=gene_names, labCol=lib_names, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
par(cex.main=1)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", Rowv=FALSE, dendrogram="column", main=main_title, ColSideColors=her2_colors, labCol=lib_names, labRow=gene_names, cexRow=0.8, cexCol=0.60, margins=c(12,6), col=rev(heat.colors(75)))
legend("bottomleft", legend=c("HER2+","HER2-"), fill=c(colors[1],colors[5]), bty="n", cex=0.8, border="white")
dev.off()


#Look at just TFAP2C and ERBB2
data_TFAP2C_ERBB2=data[TFAP2C_ERBB2_ind,]
feat_data_TFAP2C_ERBB2=feat_data[TFAP2C_ERBB2_ind,]
exp_status_TFAP2C_ERBB2=exp_status[TFAP2C_ERBB2_ind,]
rownames(data_TFAP2C_ERBB2)=feat_data_TFAP2C_ERBB2[,"Seq_Name"]
gene_names=rownames(data_TFAP2C_ERBB2)


#Heatmap - TFAP2C
#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(data_TFAP2C_ERBB2+1)
pdf(file=TFAP2C_heatmap_outfile)
x=as.matrix(z)
main_title="TFAP2C RNAseq gene-level expression"
clab=cbind(subtype_colors,her2_colors,quality_colors)
colnames(clab) = c("Subtype","Her2","Quality")
#heatmap.plus(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Rowv=NA, main=main_title, ColSideColors=clab, labRow=gene_names, labCol=lib_names, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
par(cex.main=1)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", Rowv=FALSE, dendrogram="column", main=main_title, ColSideColors=her2_colors, labCol=lib_names, labRow=gene_names, cexRow=0.8, cexCol=0.60, margins=c(12,6), col=rev(heat.colors(75)))
legend("bottomleft", legend=c("HER2+","HER2-"), fill=c(colors[1],colors[5]), bty="n", cex=0.8, border="white")
dev.off()



#Create waterfall plot
data_TFAP2C=data[TFAP2C_ind,]
data_TFAP2C_sorted=sort(data_TFAP2C, decreasing=TRUE)
her2_colors_sorted=her2_colors[order(data_TFAP2C, decreasing=TRUE)]
ymin=floor(min(data_TFAP2C))
ymax=ceiling(max(data_TFAP2C))
barplot_values=barplot(as.matrix(data_TFAP2C_sorted), plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="TFAP2C")
pdf(file=waterfall_file)
barplot(as.numeric(data_TFAP2C_sorted), ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="TFAP2C expression level", names.arg=names(data_TFAP2C_sorted), col=her2_colors_sorted)
#draw line midway between bars straddling sensitive/resistant threshold
box()
legend("topright", legend=c("ERBB2pos","ERBB2neg"), fill=unique(her2_colors_sorted))
#abline(v=barplot_values[22]+(barplot_values[23]-barplot_values[22])/2, lwd=2, col="black")
dev.off()


#Look at correlations among potential amplicon members
cor(t(data_TFAP2Camp), method="spearman")

