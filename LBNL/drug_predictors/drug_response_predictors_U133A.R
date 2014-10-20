library(randomForest)
library(mclust)
require(Hmisc)
library(ROCR)
library("gplots")

#Specify working directories and predictor data file (4 different versions) - choose one
#RMA - Standard CDF
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v1/RMA_standardCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v1/RMA_standardCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered/Neve_AffyRMA_genelevel_maxvar_stringent.csv"

#RMA - Standard CDF v2
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v2/RMA_standardCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/RMA_standardCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyRMA_genelevel_maxvar_stringent.csv"

#GCRMA - standard CDF
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v1/GCRMA_standardCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v1/GCRMA_standardCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered/Neve_AffyGCRMA_genelevel_maxvar_stringent.csv"

#GCRMA - standard CDF v2
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v2/GCRMA_standardCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/GCRMA_standardCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyGCRMA_genelevel_maxvar_stringent.csv"

#RMA - custom CDF
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v1/RMA_customCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v1/RMA_customCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered/Neve_AffyRMAcustom_genelevel_stringent.csv"

#RMA - custom CDF v2
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v2/RMA_customCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/RMA_customCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyRMAcustom_genelevel_stringent.csv"

#GCRMA custom CDF
#outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v1/GCRMA_customCDF/"
#RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v1/GCRMA_customCDF/RFmodels/"
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered/Neve_AffyGCRMAcustom_genelevel_stringent.csv"

#GCRMA custom CDF v2
outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/U133A/allDrugs/v2/GCRMA_customCDF/"
RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/GCRMA_customCDF/RFmodels/"
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyGCRMAcustom_genelevel_stringent.csv"


#Set working directory and output files
setwd(outdir)
outfile="RF_results_summary.txt"
pdf_results="RF_results.pdf"

#Import combined predictor data
raw_data_import=read.csv(datafile, row.names=1)

#Break data into features info and expression values
raw_feat_data=rownames(raw_data_import)
raw_data=raw_data_import

#Fix column names
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"

#Import cell line data
#cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list.txt"
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list_v2.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Filter down to set of cell lines with both drug data and U133A data (exclude non-malignant, non-BCCL)
#core_cell_lines=c("BT20","HCC1143","HCC1187","HCC1500","HCC1569","HCC1937","HCC1954","HCC3153","HCC70","MDAMB468","SUM149PT","BT549","HCC38","HS578T","MDAMB157","MDAMB231","MDAMB436","SUM1315MO2","SUM159PT","600MPE","AU565","BT474","BT483","CAMA1","HCC1428","HCC202","HCC2185","LY2","MCF7","MDAMB134VI","MDAMB175VII","MDAMB361","MDAMB415","MDAMB453","SKBR3","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","ZR751","ZR7530","ZR75B")
core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HBL100","HCC1143","HCC1187","HCC1428","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","ZR751","ZR7530","ZR75B")

data=raw_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#drug response data
#drugdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_10Feb.csv"
drugdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_26Aug11.csv"
raw_drugdata_import=read.csv(drugdatafile)
drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Fix row names
rownames(drug_data)[which(rownames(drug_data)=="Hs578T")]="HS578T"

#Retrieve data for only libraries in core cell line set
drug_data_filt=drug_data[core_cell_lines,]

#Transform data to -log values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#Load mean GI50 values from file
#meanGI50file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/GI50meanThresholds.csv"
meanGI50file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/GI50meanThresholds_v2.csv"
meanGI50_import=read.csv(meanGI50file, row.names=1)
drugnames=rownames(meanGI50_import)
meanGI50_data=meanGI50_import[,1]

#Fix problem drug name:
drugnames[drugnames=="Tykerb:IGF1R (1:1)"]="Tykerb(IGF1R)"

###For testing###
#drugs=drugs[1:2]
#################

#For each drug of interest, find predictors associated with response and build predictor for drug response
#First, create dataframe to hold results
predictor_perf = data.frame(cbind(drugs, N=NA, sens=NA, spec=NA, acc=NA, err=NA, res_err=NA, sen_err=NA, AUC=NA, meanGI50=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drugs

#Initialize pdf for all output
pdf(file=pdf_results)

for (i in 1:length(drugs)){ #start drug loop
drug=drugs[i]
drugname=drugnames[i]
print(paste("processing",drug))

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]

#Retrieve mean from pre-calculated file
mean_cutoff=meanGI50_data[i]

drug_data_interest_NA=which(is.na(drug_data_interest))
resistants=which(drug_data_interest<=mean_cutoff)
sensitives=which(drug_data_interest>mean_cutoff)

response_class=vector(length=length(drug_data_interest))
response_class[drug_data_interest_NA]=NA
response_class[sensitives]="sensitive"
response_class[resistants]="resistant"

#Exclude libs where response_class=NA
nonNA=which(!is.na(response_class))
data_nonNA=data[,nonNA]
response_class_nonNA=response_class[nonNA]
cell_line_data_nonNA=cell_line_data[nonNA,]

#Make sure there are at least 10 libs classified as sensitive and resistant
num_sensitive=length(which(response_class_nonNA=="sensitive"))
num_resistant=length(which(response_class_nonNA=="resistant"))

target=as.factor(response_class_nonNA)

#Down-sampling?
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

rf_model=randomForest(x=t(data_nonNA), y=target, importance = TRUE, ntree = 50001, proximity=TRUE, sampsize=sampsizes)

#Save model
rf_model_file=paste("RF_model",drugname,sep="_")
save(rf_model, file=paste(RF_model_outdir,rf_model_file,sep=""))

#Get importance measures
rf_importances=importance(rf_model, scale=FALSE)

#Determine performance statistics
confusion=rf_model$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_model$err.rate[length(rf_model$err.rate[,1]),1]*100
class1_error=paste(rownames(confusion)[1],"error rate =",format(confusion[1,3]*100, digits=4), sep=" ")
class2_error=paste(rownames(confusion)[2],"error rate =",format(confusion[2,3]*100, digits=4), sep=" ")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity =",format(sensitivity, digits=4), sep=" ")
spec_out=paste("specificity =",format(specificity, digits=4), sep=" ")
err_out=paste("overall error rate =",format(overall_error, digits=4),sep=" ")
acc_out=paste("overall accuracy =",format(overall_accuracy, digits=4),sep=" ")
misclass_1=paste(confusion[1,2], "/", num_resistant, rownames(confusion)[1],"misclassed", sep=" ")
misclass_2=paste(confusion[2,1], "/", num_sensitive, rownames(confusion)[2],"misclassed", sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Produce graph of variable importances for top 30 markers
#pdf(file=varimp_pdffile)
#varImpPlot(rf_model, type=2, n.var=30, scale=FALSE, main="Var. Imp. (Gini) for top 30 predictors")
#dev.off()

#Create ROC curve plot and calculate AUC
#Can use Sensitive vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling sensitive vs resistant
predictions=as.vector(rf_model$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
AUC_out=paste("AUC =",format(AUC, digits=5, scientific=FALSE), sep=" ")
#Then, plot the actual ROC curve
#perf_ROC=performance(pred,"tpr","fpr")
#pdf(file=ROC_pdffile)
#plot(perf_ROC, main="ROC plot")
#text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
#dev.off()

#Add performance stats to dataframe
predictor_perf[drug,"N"]=length(target)
predictor_perf[drug,"sens"]=sensitivity
predictor_perf[drug,"spec"]=specificity
predictor_perf[drug,"acc"]=overall_accuracy
predictor_perf[drug,"err"]=overall_error
predictor_perf[drug,"res_err"]=confusion[1,3]*100
predictor_perf[drug,"sen_err"]=confusion[2,3]*100
predictor_perf[drug,"AUC"]=AUC
predictor_perf[drug,"meanGI50"]=mean_cutoff

#Create side by side plot with MDS (or vote distribution) plot and performance stats
#Produce MDS plot
#save current par settings
op = par(no.readonly = TRUE)
par(mfrow=c(1,2), oma=c(4,2,10,2)) #bottom,left,top,right
par(mar=c(10,4,0,0))
#target_labels=as.vector(target)
#target_labels[target_labels=="sensitive"]="S"
#target_labels[target_labels=="resistant"]="R"
#MDSplot(rf_model, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")

#Produce back-to-back histogram of vote distributions for Sensitive and Resistant
options(digits=2) 
out <- histbackback(split(rf_model$votes[,"sensitive"], target), probability=FALSE, xlim=c(-15,15), axes=TRUE, ylab="Fraction votes (sensitive)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 

#Add stats
par(mar=c(4,1,0,1))
plot.new()
stats_legend=c(
"Classification Stats:",
paste("N =",length(target), "cell lines", sep=" "),
sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2,AUC_out
)
legend("topleft", legend=stats_legend, bty="n")
title(main=drugname, outer=TRUE, line=2)

#Set back to old par settings
par(op)

#Create heatmap of top X predictors
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

top30=sort(rf_importances[,4], decreasing=TRUE)[1:30]
top30_data=data_nonNA[names(top30),cell_lines_sorted]
predictor_names=rownames(top30_data)
lib_names=colnames(top30_data)

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[cell_lines_sorted,"BCCLclassification2"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="ERBB2Amp"]="blue"
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"
subtype_colors[subtype_colors=="Non-malignant"]="yellow"
subtype_colors[subtype_colors=="Unknown"]="pink"

main_title=drugname
heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, labCol=lib_names, labRow=predictor_names, colsep=num_sensitive, sepcolor="black", cexRow=0.7, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))
#heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, labCol=lib_names, labRow=predictor_names, colsep=num_sensitive, sepcolor="black", cexRow=0.7, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))
legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)

#Create corresponding histogram/waterfall plot for drug response data
ymin=floor(min(drug_data_interest_sorted))
ymax=ceiling(max(drug_data_interest_sorted))
barplot_values=barplot(drug_data_interest_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drugname)
barplot(drug_data_interest_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
#draw line midway between bars straddling sensitive/resistant threshold
cutoff_line=barplot_values[num_sensitive]+(barplot_values[num_sensitive+1]-barplot_values[num_sensitive])/2
abline(v=cutoff_line, lwd=2, col="black")
legend("topright", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), cex=0.9)
text(x=cutoff_line-(length(target)*0.3), y=ymax-((ymax-ymin)*0.05), labels=paste("cutoff =", format(mean_cutoff, digits=4), sep=" "), pos=4)

} #end drug loop
dev.off()

#Write summary table to file
write.table(predictor_perf, file=outfile, sep="\t", row.names=FALSE)

