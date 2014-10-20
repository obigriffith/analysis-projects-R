library("gplots")
library("heatmap.plus")
library(randomForest)
library(ROCR)

#Cell line info
#celllinefile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/cell_line_info.txt"
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
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/analysis/ERBB2/"
#outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/analysis/ERBB2/RNA_vs_exonArray/"

setwd(outdir)

#HER2amp_heatmap_outfile = "BCCL_HER2Amp_67libs_heatmap.pdf"
HER2amp_heatmap_outfile = "BCCL_HER2Amp_67libs_heatmap_extended.pdf"
outfile="RFoutput.txt"
case_pred_outfile="CasePredictions.txt"
varimp_pdffile="varImps.pdf"
MDS_pdffile="MDS.pdf"
case_margins_file="Margins.pdf"
ROC_pdffile="ROC.pdf"
vote_dist_pdffile="vote_dist.pdf"
waterfallfile="waterfall.pdf"

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
#ERBB2=as.vector(cell_data[,"ERBB2"])
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

#HER2 status
her2_colors=ERBB2
her2_colors[her2_colors=="HighAmp"]=colors[1]
her2_colors[her2_colors=="Amp"]=colors[1]
her2_colors[her2_colors=="OverExp"]=colors[2]
her2_colors[her2_colors=="LowAmp"]=colors[3]
her2_colors[her2_colors=="NoAmp"]=colors[5]
her2_colors[her2_colors=="ND"]=colors[6]

#Quality status
quality_colors=Qualities
quality_colors[quality_colors=="High"]=colors[1]
quality_colors[quality_colors=="Int"]=colors[3]
quality_colors[quality_colors=="Low"]=colors[5]

#Get ERBB2 Amp genes
### exp set - flanking 1M bp (500,000bp up- and down-stream ###
FBXL20_ind=which(feat_data[,"Seq_Name"]=="FBXL20")
PPARBP_ind=which(feat_data[,"Seq_Name"]=="PPARBP") #aka MED1
CRKRS_ind=which(feat_data[,"Seq_Name"]=="CRKRS")
### exp set ### - core amplicon
NEUROD2_ind=which(feat_data[,"Seq_Name"]=="NEUROD2")
PPP1R1B_ind=which(feat_data[,"Seq_Name"]=="PPP1R1B")
STARD3_ind=which(feat_data[,"Seq_Name"]=="STARD3")
TCAP_ind=which(feat_data[,"Seq_Name"]=="TCAP")
PNMT_ind=which(feat_data[,"Seq_Name"]=="PNMT")
PERLD1_ind=which(feat_data[,"Seq_Name"]=="PERLD1")
ERBB2_ind=which(feat_data[,"Seq_Name"]=="ERBB2")
C17orf37_ind=which(feat_data[,"Seq_Name"]=="C17orf37")
GRB7_ind=which(feat_data[,"Seq_Name"]=="GRB7")
IKZF3_ind=which(feat_data[,"Seq_Name"]=="IKZF3")
### exp set ### - core amplicon
ZPBP2_ind=which(feat_data[,"Seq_Name"]=="ZPBP2")
GSDMB_ind=which(feat_data[,"Seq_Name"]=="GSDMB")
ORMDL3_ind=which(feat_data[,"Seq_Name"]=="ORMDL3")
GSDMA_ind=which(feat_data[,"Seq_Name"]=="GSDMA")
PSMD3_ind=which(feat_data[,"Seq_Name"]=="PSMD3")
THRAP4_ind=which(feat_data[,"Seq_Name"]=="THRAP4") #aka DRIP100, MED24
THRA_ind=which(feat_data[,"Seq_Name"]=="THRA")
NR1D1_ind=which(feat_data[,"Seq_Name"]=="NR1D1")
MSL1_ind=which(feat_data[,"Seq_Name"]=="MSL1")
CASC3_ind=which(feat_data[,"Seq_Name"]=="CASC3")
RAPGEFL1_ind=which(feat_data[,"Seq_Name"]=="RAPGEFL1")
### exp set ### - flanking 1M bp (500,000bp up- and down-stream ###

#AmpGenes_ind=c(NEUROD2_ind,PPP1R1B_ind,STARD3_ind,TCAP_ind,PNMT_ind,PERLD1_ind,ERBB2_ind,C17orf37_ind,GRB7_ind,IKZF3_ind)
AmpGenes_ind=c(FBXL20_ind,PPARBP_ind,CRKRS_ind,NEUROD2_ind,PPP1R1B_ind,STARD3_ind,TCAP_ind,PNMT_ind,PERLD1_ind,ERBB2_ind,C17orf37_ind,GRB7_ind,IKZF3_ind,ZPBP2_ind,GSDMB_ind,ORMDL3_ind,GSDMA_ind,PSMD3_ind,THRAP4_ind,THRA_ind,NR1D1_ind,MSL1_ind,CASC3_ind,RAPGEFL1_ind)
AmpGenes_ind2=c(STARD3_ind,PERLD1_ind,ERBB2_ind,C17orf37_ind,GRB7_ind)

data_ERBB2amp=data[AmpGenes_ind,]
feat_data_ERBB2amp=feat_data[AmpGenes_ind,]
exp_status_ERBB2amp=exp_status[AmpGenes_ind,]
rownames(data_ERBB2amp)=feat_data_ERBB2amp[,"Seq_Name"]
gene_names=rownames(data_ERBB2amp)

#Calculate ranks, use just 5 "most useful" genes in amplicon?
data_ERBB2amp2=data[AmpGenes_ind2,]
data_ERBB2amp2_ranks=t(apply(data_ERBB2amp2,1,rank))
data_ERBB2amp2_ranksums=apply(data_ERBB2amp2_ranks,2,sum)
cell_lines_sorted=names(sort(data_ERBB2amp2_ranksums, decreasing=TRUE)) 
data_ERBB2amp2_ranksums_sorted=data_ERBB2amp2_ranksums[cell_lines_sorted]
subtype_colors_sorted=subtype_colors[cell_lines_sorted]

#Create waterfall plot using rank sums
ymin=floor(min(data_ERBB2amp2_ranksums_sorted))
ymax=ceiling(max(data_ERBB2amp2_ranksums_sorted))
barplot_values=barplot(data_ERBB2amp2_ranksums_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="HER2amp ranksum")
pdf(file=waterfallfile)
barplot(data_ERBB2amp2_ranksums_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="HER2amp ranksum", col=subtype_colors_sorted)
#draw line midway between bars straddling sensitive/resistant threshold
box()
legend("topright", legend=c("Luminal","Unknown","Basal","Claudin-low","Normal-like","Normal"), fill=unique(subtype_colors_sorted))
abline(v=barplot_values[22]+(barplot_values[23]-barplot_values[22])/2, lwd=2, col="black")
dev.off()







#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(data_ERBB2amp+1)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap
pdf(file=HER2amp_heatmap_outfile)
x=as.matrix(z)
main_title="HER2 Amplicon RNAseq gene-level expression"
#main_title="RNAseq transcript-level expression"
clab=cbind(subtype_colors,her2_colors,quality_colors)
colnames(clab) = c("Subtype","Her2","Quality")
#heatmap.plus(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Rowv=NA, main=main_title, ColSideColors=clab, labRow=gene_names, labCol=lib_names, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
par(cex.main=1)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", Rowv=FALSE, dendrogram="column", main=main_title, ColSideColors=her2_colors, labCol=lib_names, labRow=gene_names, cexRow=0.8, cexCol=0.60, margins=c(12,6), col=rev(heat.colors(75)))
#legend("bottomleft", legend=c("HighAmp","OverExp","LowAmp","NoAmp","ND"), fill=c(colors[1],colors[2],colors[3],colors[5],colors[6]), bty="n", cex=0.8, border="white")
legend("bottomleft", legend=c("Amp","NoAmp"), fill=c(colors[1],colors[5]), bty="n", cex=0.8, border="white")
dev.off()

#RF analysis
#Make binary Amp/NoAmp ERBB2 factor
ERBB2_binary=ERBB2
ERBB2_binary[ERBB2_binary=="HighAmp"]="Amp"
ERBB2_binary[ERBB2_binary=="LowAmp"]="NoAmp"
ERBB2_binary[ERBB2_binary=="ND"]="NoAmp"
ERBB2_binary[ERBB2_binary=="OverExp"]="NoAmp"
target=ERBB2_binary

#attempt to classify by subtype
rf_output=randomForest(x=t(z), y=as.factor(target), importance = TRUE, ntree = 50001, proximity=TRUE)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=100-overall_error

#Prepare stats for output to file
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(acc_out,err_out), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
target_labels=as.vector(target)
target_labels[target_labels=="Amp"]="A"
target_labels[target_labels=="NoAmp"]="N"

pdf(file=MDS_pdffile)
MDSplot(rf_output, as.factor(target_labels), k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
dev.off()

#Create ROC curve plot and calculate AUC
#Can use vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling one class vs other
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()

