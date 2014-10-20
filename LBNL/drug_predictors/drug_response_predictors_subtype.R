library(randomForest)
library(mclust)
require(Hmisc)
library(ROCR)
library("gplots")

outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/subtype/"
RF_model_outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/subtype/RFmodels/"
setwd(outdir)
outfile="RF_results_summary.txt"
pdf_results="RF_results.pdf"

#Import cell line data
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list_v2.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Import subtype data
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/SubtypeVariables_v2.txt"
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)

#Exclude lines with subtype of normal-like or unknown
raw_data=raw_data_import[which(raw_data_import[,"Luminal"]==1 | raw_data_import[,"Basal"]==1 | raw_data_import[,"Claudin.low"]==1),]

#drug response data
drugdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_26Aug11.csv"
raw_drugdata_import=read.csv(drugdatafile, na.strings = c("NaN","","NA","N/A"))
drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Fix row names
rownames(drug_data)[which(rownames(drug_data)=="Hs578T")]="HS578T"

#Filter down to core set of cell lines - must have drug data and at least one other molecular profiling data type
core_cell_lines=rownames(raw_data)[which(rownames(raw_data) %in% rownames(drug_data))]
data=raw_data[core_cell_lines,]
cell_line_data=raw_cell_line_import[core_cell_lines,]
drug_data_filt=drug_data[core_cell_lines,]

#Transform data to -log10 values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#Retrieve mean GI50 values
GI50meanfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/GI50meanThresholds_v2.csv"
GI50mean_import=read.csv(GI50meanfile)
GI50mean_data=GI50mean_import[,2]
names(GI50mean_data)=GI50mean_import[,1]
drug_names=as.vector(GI50mean_import[,1])
cbind(drugs,names(GI50mean_data)) #Check to make sure drugs are in same order in drug data and mean GI50 files

#Fix problem drug name:
drug_names[drug_names=="Tykerb:IGF1R (1:1)"]="Tykerb(IGF1R)"

#For each drug of interest, find predictors associated with response and build predictor for drug response
#First, create dataframe to hold results
predictor_perf = data.frame(cbind(drugs, N=NA, sens=NA, spec=NA, acc=NA, err=NA, res_err=NA, sen_err=NA, AUC=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drugs

#Initialize pdf for all output
pdf(file=pdf_results)

for (i in 1:length(drugs)){ #start drug loop
#for (i in 1:2){ #TESTING

drug=drugs[i]
drug_name=drug_names[i]
print(paste("processing",drug_name))

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]

mean_cutoff=GI50mean_data[i]
drug_data_interest_NA=which(is.na(drug_data_interest))
resistants=which(drug_data_interest<=mean_cutoff)
sensitives=which(drug_data_interest>mean_cutoff)

response_class=vector(length=length(drug_data_interest))
response_class[drug_data_interest_NA]=NA
response_class[sensitives]="sensitive"
response_class[resistants]="resistant"

#Exclude libs where response_class=NA
nonNA=which(!is.na(response_class))
data_nonNA=data[nonNA,]
response_class_nonNA=response_class[nonNA]
cell_line_data_nonNA=cell_line_data[nonNA,]
cell_lines_nonNA=rownames(data_nonNA)

#Make sure there are at least 10 libs classified as sensitive and resistant
num_sensitive=length(which(response_class_nonNA=="sensitive"))
num_resistant=length(which(response_class_nonNA=="resistant"))

target=as.factor(response_class_nonNA)

#Down-sampling
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

####MAKE SURE FACTORS TREATED AS FACTORS####
data_nonNA=data.frame(lapply(data_nonNA, as.factor))
rownames(data_nonNA)=cell_lines_nonNA

#Run RF on imputed data. Recall, you must extract all but first column where rfImpute inserts target there
#TESTING
#rf_model=randomForest(x=data_imputed[,2:length(colnames(data_imputed))], y=target, importance = TRUE, ntree = 501, proximity=TRUE, sampsize=sampsizes)
#rf_model=randomForest(x=data_imputed[,2:length(colnames(data_imputed))], y=target, importance = TRUE, ntree = 10001, proximity=TRUE)
#rf_model=randomForest(x=data_imputed[,2:length(colnames(data_imputed))], y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)
rf_model=randomForest(x=data_nonNA, y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)

#Save model
rf_model_file=paste("RF_model",drug_name,sep="_")
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


#Plot error rates as function of number of trees
plot(rf_model, main="")
legend("topright", legend=err_out, bty="n")

#Produce graph of variable importances for top 30 markers
#pdf(file=varimp_pdffile)
#varImpPlot(rf_model, type=2, n.var=30, scale=FALSE, main="Var. Imp. (Gini) for top 30 predictors")
#dev.off()

#Produce back-to-back histogram of vote distributions for Sensitive and Resistant
#options(digits=2) 
#pdf(file=vote_dist_pdffile)
#out <- histbackback(split(rf_model$votes[,"sensitive"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for cell lines classified by RF', axes=TRUE, ylab="Fraction votes (sensitive)")
#add color
#barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
#barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
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

#Create side by side plot with MDS plot and performance stats
#Produce MDS plot
#save current par settings
op = par(no.readonly = TRUE)

par(mfrow=c(1,2), mar=c(0.5,0.5,2,0), oma=c(4,4,4,4))
target_labels=as.vector(target)
target_labels[target_labels=="sensitive"]="S"
target_labels[target_labels=="resistant"]="R"
MDSplot(rf_model, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")

#Add stats
plot.new()
stats_legend=c(
"Classification Stats",
paste("n =",length(target), "cell lines", sep=" "),
sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2,AUC_out
)
legend("left", legend=stats_legend, bty="n")
title(main=drug, outer=TRUE)

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

topX_data=data_nonNA[cell_lines_sorted,]
predictor_names=colnames(topX_data)
lib_names=rownames(topX_data)

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[cell_lines_sorted,"BCCLclassification2"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="ERBB2Amp"]="blue"
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"
subtype_colors[subtype_colors=="Non-malignant"]="yellow"
subtype_colors[subtype_colors=="Unknown"]="pink"

main_title=drug
#pdf(file=heatmapfile)
#heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=CleanNames, colsep=num_sensitive, sepcolor="black", cexRow=0.8, cexCol=0.75, margins=c(6,12), col=rev(heat.colors(75)))
#legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","ExomeSeq","ManMut","CNV","RPPA","WB","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.7)
#legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)
#dev.off()

z=data.frame(lapply(topX_data, as.numeric))-1

heatmap.2(as.matrix(t(z)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, symbreaks=FALSE, key=FALSE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, labCol=lib_names, labRow=predictor_names, colsep=num_sensitive, sepcolor="black", cexRow=0.7, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))

legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)
#dev.off()

#Create corresponding histogram/waterfall plot for drug response data
ymin=floor(min(drug_data_interest_sorted))
ymax=ceiling(max(drug_data_interest_sorted))
barplot_values=barplot(drug_data_interest_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
#pdf(file=waterfallfile)
barplot(drug_data_interest_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
#draw line midway between bars straddling sensitive/resistant threshold
cutoff_line=barplot_values[num_sensitive]+(barplot_values[num_sensitive+1]-barplot_values[num_sensitive])/2
abline(v=cutoff_line, lwd=2, col="black")
legend("topright", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), cex=0.9)
text(x=cutoff_line+1, y=ymax-1, labels=paste("cutoff =", format(mean_cutoff, digits=4), sep=" "), pos=4)
#dev.off()

} #end drug loop
dev.off()

#Write summary table to file
write.table(predictor_perf, file=outfile, sep="\t", row.names=FALSE)

