library(randomForest)
library("gplots")

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/GreenMouseModels/"
setwd(outdir)
heatmapfile="breastExon_genelevel_stringent_50lines_unsup_heatmap.pdf"
MDS_pdffile="breastExon_genelevel_stringent_50lines_MDS.pdf"
outfile="breastExon_genelevel_stringent_50lines_RFresults.txt"

#Import cell line data
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

### exonarray data ###
EAdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ExonArray/filtered/breastExon_genelevel_stringent.csv"
EA_data_import=read.csv(EAdatafile) 
EA_feat_data=as.vector(EA_data_import[,1])
EA_data=EA_data_import[,2:length(colnames(EA_data_import))]
rownames(EA_data)=EA_feat_data

#Fix misnamed libraries
colnames(EA_data)[which(colnames(EA_data)=="X184A1")]="184A1"
colnames(EA_data)[which(colnames(EA_data)=="X184B5")]="184B5"
colnames(EA_data)[which(colnames(EA_data)=="X600MPE")]="600MPE"

#Define set of cell lines for analysis - all Basal, Luminal, Claudin-low (excluding "outlier", Normal-like, "unknown")
core_cell_lines=c("BT20","HCC1143","HCC1187","HCC1500","HCC1599","HCC1806","HCC1937","HCC1954","HCC3153","HCC70","MDAMB468","SUM102PT","SUM149PT","BT549","HCC1395","HCC38","HS578T","MDAMB157","MDAMB231","MDAMB436","SUM1315MO2","SUM159PT","600MPE","AU565","BT474","BT483","CAMA1","HCC1419","HCC1428","HCC202","HCC2185","HCC2218","LY2","MCF7","MDAMB134VI","MDAMB175VII","MDAMB361","MDAMB415","MDAMB453","SKBR3","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")
data=EA_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#Create heatmap to display unsupervised clustering of lines
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[,"BCCLclassification"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"

pdf(file=heatmapfile)
main_title="Exon-array, stringent filtered, unsupervised clustering"
heatmap.2(as.matrix(data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="column", ColSideColors=subtype_colors, labRow=FALSE, cexRow=0.7, cexCol=0.9, cex.main=0.7, margins=c(6,4), col=rev(heat.colors(75)))
legend("bottomleft", legend=c("Luminal","Basal","Claudin-low"), fill=c("black","red","green"), title="Subtype", border="white",  bty="n",cex=0.7)
legend("left",legend=c("genes=1364","lines=50"), border="white", bty="n", cex=0.7)
dev.off()


#Use Randomforests to classify according to subtype
target=as.factor(as.vector(cell_line_data[,"BCCLclassification"]))

#Transpose data so that predictor variables are columns
tdata=as.data.frame(t(data))

#Run RF
rf_model=randomForest(x=tdata, y=target, importance = TRUE, ntree = 50001, proximity=TRUE)

#Save model
rf_model_file="breastExon_genelevel_stringent_50lines_RFmodel"
save(rf_model, file=rf_model_file)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
class1_error=paste(rownames(confusion)[1],"error rate =",format(confusion[1,4]*100, digits=4), sep=" ")
class2_error=paste(rownames(confusion)[2],"error rate =",format(confusion[2,4]*100, digits=4), sep=" ")
class3_error=paste(rownames(confusion)[3],"error rate =",format(confusion[3,4]*100, digits=4), sep=" ")
overall_accuracy=100-overall_error

#Prepare stats for output to file
err_out=paste("overall error rate =",format(overall_error, digits=4),sep=" ")
acc_out=paste("overall accuracy =",format(overall_accuracy, digits=4),sep=" ")
misclass_1=paste(confusion[1,2]+confusion[1,3], "/", sum(confusion[1,1:3]), rownames(confusion)[1],"misclassed", sep=" ")
misclass_2=paste(confusion[2,1]+confusion[2,3], "/", sum(confusion[2,1:3]), rownames(confusion)[2],"misclassed", sep=" ")
misclass_3=paste(confusion[3,1]+confusion[3,2], "/", sum(confusion[3,1:3]), rownames(confusion)[3],"misclassed", sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:3,1:3]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Create side by side plot with MDS plot and performance stats
#Produce MDS plot
pdf(file=MDS_pdffile)
par(mfrow=c(1,2), mar=c(0.5,0.5,2,0), oma=c(4,4,4,4))
target_labels=as.vector(target)
target_labels[target_labels=="Basal"]="B"
target_labels[target_labels=="Claudin-low"]="C"
target_labels[target_labels=="Luminal"]="L"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "green", "black"), main="MDS plot")

#Add stats
plot.new()
stats_legend=c(
"Classification Stats",
paste("n =",length(target), "cell lines", sep=" "),
acc_out,err_out,class1_error,class2_error,class3_error,misclass_1,misclass_2,misclass_3)
legend("left", legend=stats_legend, bty="n")
dev.off()

#Print results to file
write.table(sort(rf_importances[,4], decreasing=TRUE),file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(acc_out,err_out,class1_error,class2_error,class3_error,misclass_1,misclass_2,misclass_3), file=outfile, append=TRUE)


#Recreate RF model with just top genes
topGenes=names(which(rf_importances[,4]>0.001))

tdata_topgenes=tdata[,topGenes]

#Run RF
rf_model2=randomForest(x=tdata_topgenes, y=target, importance = TRUE, ntree = 50001, proximity=TRUE)

#Save model
rf_model_file2="breastExon_genelevel_top115genes_50lines_RFmodel"
save(rf_model2, file=rf_model_file2)

