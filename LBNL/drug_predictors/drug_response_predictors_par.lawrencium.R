#!/usr/bin/env Rscript

args=(commandArgs(TRUE))
drug_i=as.numeric(args[1])
outdir_base=args[2]
datatypes_string=args[3]
datatypes=as.vector(strsplit(datatypes_string,",")[[1]])
datafile=args[4]
core_cell_lines_string=args[5]
core_cell_lines=as.vector(strsplit(core_cell_lines_string,",")[[1]])

#Command line usage:
#/global/home/users/olgriff/ac_cancerseq/drug_predictors/drug_response_predictors_par.R 1 /global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/mutation/ "MM" "600MPE,AU565,BT20,BT474,BT483,BT549,CAMA1,HCC1143,HCC1187,HCC1395,HCC1419,HCC1569,HCC1806,HCC1937,HCC1954,HCC38,HCC70,HS578T,MCF7,MDAMB134VI,MDAMB157,MDAMB175VII,MDAMB231,MDAMB361,MDAMB415,MDAMB436,MDAMB453,MDAMB468,SKBR3,SUM1315MO2,SUM149PT,SUM159PT,SUM185PE,SUM52PE,T47D,UACC812,UACC893,ZR751,ZR7530"

#print(drug_i)

library(randomForest)
#library(mclust)
#require(Hmisc)
library(ROCR)
library("gplots")

#data files - HARD CODED
#cell_line_datafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/BCCL_data_list.txt"
cell_line_datafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/BCCL_data_list_v2.txt"
#datafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/combined/BCCL_combined_data.3.COSMICextended.imp.txt"
#datafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/combined/BCCL_combined_data.4.imp.txt"
#datafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/combined/BCCL_combined_data.4.txt"
#drugdatafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_10Feb.csv"
drugdatafile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_26Aug11.csv"
GI50meanfile="/global/home/users/olgriff/ac_cancerseq/drug_predictors/drugdata/GI50meanThresholds_v2.csv"

#output folders/files
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/core_alldata_imputed/"
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/rnaseq/" #RNAseq - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/exonarray/" #exon array - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/U133A/" #U133A - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/Meth/" #Meth - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/mutation/" #mutation - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/RPPA/" #RPPA - all lines
#outdir_base="/global/home/users/olgriff/ac_cancerseq/drug_predictors/RFpredictors/allDrugs/independent/SNP/" #SNP6 - all lines

RF_model_outdir=paste(outdir_base,"RFmodels/",sep="")
RF_pdf_outdir=paste(outdir_base,"pdfs/",sep="")
RF_results_outdir=paste(outdir_base,"results/",sep="")
setwd(outdir_base)

#Set datatypes to include in analysis
#datatypes=c("RS","EA","U133A","Meth","MM","RPPA","SNP") #all
#datatypes=c("RS") #RNAseq
#datatypes=c("EA") #Exon Array
#datatypes=c("U133A") #U133A
#datatypes=c("Meth") #Meth
#datatypes=c("MM") #Mutations
#datatypes=c("RPPA") #RPPA
#datatypes=c("SNP") #SNP6

#Set core cell lines - this can change depending on which input data and which analysis
#Choose One
#core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")

#RNAseq lines
#core_cell_lines=c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE","AU565","BT474","BT483","BT549","CAMA1","EFM192A","EFM192B","EFM192C","HCC1143","HCC1395","HCC1419","HCC1428","HCC1569","HCC1599","HCC1806","HCC1937","HCC1954","HCC202","HCC2218","HCC3153","HCC38","HCC70","HS578T","JIMT1","LY2","MB157","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","MX1","SKBR3","SUM1315MO2","SUM149PT","SUM225CWN","SUM229PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")

#Exon-array lines
#core_cell_lines=c("184A1","184B5","600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HBL100","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1569","HCC1599","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC2218","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM102PT","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")

#U133A lines
#core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HBL100","HCC1143","HCC1187","HCC1428","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","ZR751","ZR7530","ZR75B")

#Meth lines
#core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HBL100","HCC1143","HCC1187","HCC1428","HCC1569","HCC1599","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","ZR751","ZR7530","ZR75B")

#Mutation lines
#core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1569","HCC1806","HCC1937","HCC1954","HCC38","HCC70","HS578T","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530")

#RPPA lines
#core_cell_lines=c("184A1","184B5","600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HBL100","HCC1143","HCC1395","HCC1419","HCC1428","HCC1569","HCC1599","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC2218","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB157","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM102PT","SUM1315MO2","SUM149PT","SUM159PT","SUM225CWN","T47D","UACC812","UACC893","ZR7530","ZR75B")

#SNP lines
#core_cell_lines=c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE","AU565","BT20","BT474","BT483","BT549","CAL120","CAL148","CAL51","CAL851","CAMA1","COLO824","EFM19","EFM192A","EFM192B","EFM192C","EVSAT","HBL100","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1569","HCC1599","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC2218","HCC3153","HCC38","HCC70","HDQP1","HS578T","JIMT1","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB157","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","MFM223","MT3","PMC42","S1","SKBR3","SUM102PT","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM229PE","SUM44PE","SUM52PE","T4","T47D","UACC812","ZR751","ZR75B")


#Retrieve mean GI50 values and drug names
GI50mean_import=read.csv(GI50meanfile)
GI50mean_data=GI50mean_import[,2]
names(GI50mean_data)=GI50mean_import[,1]
drug_names=as.vector(GI50mean_import[,1])

#drug response data
raw_drugdata_import=read.csv(drugdatafile, na.strings = c("NaN","","NA","N/A"))
drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Fix row names
rownames(drug_data)[which(rownames(drug_data)=="Hs578T")]="HS578T"

###Limit cell lines to overlap with drug data
core_cell_lines=core_cell_lines[which(core_cell_lines %in% rownames(drug_data))]

#Retrieve data for only libraries in core cell line set
drug_data_filt=drug_data[core_cell_lines,]

#Transform drug data to -log10 values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

cbind(drugs,names(GI50mean_data)) #Check to make sure drugs are in same order in drug data and mean GI50 files

#For drug of interest, find predictors associated with response and build predictor for drug response
drug=drugs[drug_i]
drug_name=drug_names[drug_i]
print(paste("processing",drug_name))

#Create output files
#dataframe to hold results
predictor_perf = data.frame(cbind(drug_name, N=NA, sens=NA, spec=NA, acc=NA, err=NA, res_err=NA, sen_err=NA, AUC=NA, opt_feat_num=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drug

#RF results txt file - performance stats
results_file=paste(drug_name,".txt",sep="")
results_path=paste(RF_results_outdir,results_file,sep="")

#RF pdf file - figures
pdf_file=paste(drug_name,".pdf",sep="")
pdf_path=paste(RF_pdf_outdir,pdf_file,sep="")

#RF model files/paths
rf_model_file=paste("RF_model",drug_name,sep="_")
rf_model_path=paste(RF_model_outdir,rf_model_file,sep="")
rf_model_file_opt=paste("RF_model_opt",drug_name,sep="_")
rf_model_opt_path=paste(RF_model_outdir,rf_model_file_opt,sep="")

#Initialize pdf for all output
pdf(file=pdf_path)

#Import cell line data
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Import predictor data
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:2]
raw_data=raw_data_import[,3:length(colnames(raw_data_import))]

#Fix column names
colnames(raw_data)[which(colnames(raw_data)=="X184A1")]="184A1"
colnames(raw_data)[which(colnames(raw_data)=="X184B5")]="184B5"
colnames(raw_data)[which(colnames(raw_data)=="X21MT1")]="21MT1"
colnames(raw_data)[which(colnames(raw_data)=="X21MT2")]="21MT2"
colnames(raw_data)[which(colnames(raw_data)=="X21NT")]="21NT"
colnames(raw_data)[which(colnames(raw_data)=="X21PT")]="21PT"
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"

#Filter down to core set of cell lines 
data=raw_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#Set row names as datatype__featureID
rownames(data)=paste(raw_feat_data[,"DataType"],raw_feat_data[,"ID"], sep="__")

#Filter down to just specified datatypes
data=data[which(raw_feat_data[,1] %in% datatypes),]
feat_data=raw_feat_data[which(raw_feat_data[,1] %in% datatypes),]
feat_IDs=feat_data[,"ID"]

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]

mean_cutoff=GI50mean_data[drug_name]
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

#Impute missing values. So far, only cell lines where drug data is NA have been excluded. Predictor data still contains NAs
#data_imputed=rfImpute(x=t(data_nonNA), y=target, ntree=300, iter=5)#Increase numbers of ntree and iter for final result?
#data_imputed=rfImpute(x=t(data_nonNA), y=target, ntree=50, iter=2)#TESTING

#Down-sampling
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

#Run RF
#rf_model=randomForest(x=t(data_nonNA), y=target, importance = TRUE, ntree = 501, proximity=TRUE, sampsize=sampsizes) #TESTING
rf_model=randomForest(x=t(data_nonNA), y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)

#Save model
save(rf_model, file=rf_model_path)

#Show performance of models with sequentially reduced numbers of variables
#rfcv_result=rfcv(t(data_nonNA), target, cv.fold=5, scale="log", step=0.9, mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE)
#with(rfcv_result, plot(n.var, error.cv, log="x", type="o", lwd=2))
#result <- replicate(3, rfcv(t(data_nonNA), target, cv.fold=5, scale="log", step=0.1, mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE), simplify=FALSE) #TESTING
result <- replicate(50, rfcv(t(data_nonNA), target, cv.fold=5, scale="log", step=0.9, mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
mean.error.cv=signif(rowMeans(error.cv), digits=2) #Round errors to two sig digits. If different features numbers both result in 30% error, then choose smaller
min_error=min(mean.error.cv[-length(mean.error.cv)]) #Find the smallest error, excluding the error for num_features==1
opt_feat_num=min(as.numeric(names(which(mean.error.cv[-length(mean.error.cv)]==min_error)))) #Take the smallest feature number with smallest error, again excluding the error for num_features=1

#Future method to try, optimize mtry?
#tuneRF(t(data_nonNA), target, ntreeTry=1000, stepFactor=1.5, improve=0, trace=TRUE, plot=TRUE, doBest=FALSE)

#Get importance measures
rf_importances=importance(rf_model, scale=FALSE)

#create and save model with optimum feature number?
opt_features=names(sort(rf_importances[,"MeanDecreaseGini"], decreasing=TRUE)[1:opt_feat_num])
opt_feat_data_nonNA=data_nonNA[opt_features,]
rf_model_opt=randomForest(x=t(opt_feat_data_nonNA), y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)
save(rf_model_opt, file=rf_model_opt_path)

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

#save current par settings, and create temp new ones
op = par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(4,4,2,1), oma=c(6,1,6,1))

#Plot error rates as function of number of trees
plot(rf_model, main="")
legend("topright", legend=err_out, bty="n")

#Plot error rates as function of numbers of features
matplot(result[[1]]$n.var, cbind(error.cv, rowMeans(error.cv)), type="l",lwd=c(rep(1, ncol(error.cv)),2), col=c(rep("black", dim(error.cv)[2]),"red"), lty=1, log="x",xlab="number of features", ylab="CV Error")
legend("topright", legend=c(paste("opt_feat_num =", opt_feat_num),paste("min_error =",format(min_error, digits=4))), bty="n")
title(main=drug, outer=TRUE)

#Set back to old par settings
par(op)

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
predictor_perf[drug,"opt_feat_num"]=opt_feat_num

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

#For heatmap purposes, use either the opt_feat_num or 1000, whichever is smaller 
#This is necessary because occasionally opt_feat_num will be huge and cause hclust to crash (insufficient memory)
if (opt_feat_num<1000){
 topX=sort(rf_importances[,4], decreasing=TRUE)[1:opt_feat_num]
}else{
 topX=sort(rf_importances[,4], decreasing=TRUE)[1:1000]
}

topX_data=data_nonNA[names(topX),cell_lines_sorted]
predictor_names=rownames(topX_data)
lib_names=colnames(topX_data)

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[cell_lines_sorted,"BCCLclassification2"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="ERBB2Amp"]="blue"
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"
subtype_colors[subtype_colors=="Non-malignant"]="yellow"
subtype_colors[subtype_colors=="Unknown"]="pink"

#Get data_type from predictor names
getPrefix=function(x){y=x[[1]][1]; return(y)}
predictor_names_split=strsplit(x=predictor_names, split="__")
prefixes=sapply(predictor_names_split,getPrefix, simplify = TRUE)

#Make predictor names without prefix, drop leading ID for RNAseq Names (only type to have 3 components)
getName=function(x){y=x[2:length(x)]; return(y)}
Names=sapply(predictor_names_split,getName,simplify=FALSE)
pasteNames=function(x){
 if (length(x)==3){
  y=paste(x[2:3],collapse="__");
 }else{
  y=paste(x,collapse="__");
 }
return(y)
}
CleanNames=sapply(Names,pasteNames)

#Set colors for datatype sidebar
colors=rainbow(7)
datatype_colors=prefixes
datatype_colors[datatype_colors=="RS"]=colors[1]
datatype_colors[datatype_colors=="EA"]=colors[2]
datatype_colors[datatype_colors=="U133A"]=colors[3]
#datatype_colors[datatype_colors=="ES"]=colors[4]
datatype_colors[datatype_colors=="MM"]=colors[4]
#datatype_colors[datatype_colors=="CNV"]=colors[6]
datatype_colors[datatype_colors=="SNP"]=colors[5]
datatype_colors[datatype_colors=="RPPA"]=colors[6]
#datatype_colors[datatype_colors=="WB"]=colors[8]
datatype_colors[datatype_colors=="Meth"]=colors[7]

main_title=drug
#pdf(file=heatmapfile)
#heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=CleanNames, colsep=num_sensitive, sepcolor="black", cexRow=0.8, cexCol=0.75, margins=c(6,12), col=rev(heat.colors(75)))
#legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","ExomeSeq","ManMut","CNV","RPPA","WB","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.7)
#legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)
#dev.off()

#pdf(file=heatmapfile_scaled)
if (length(CleanNames)<=50){
 heatmap.2(as.matrix(topX_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=CleanNames, colsep=num_sensitive, sepcolor="black", cexRow=0.6, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))
}else{
 heatmap.2(as.matrix(topX_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=FALSE, colsep=num_sensitive, sepcolor="black", cexRow=0.7, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))
}

#legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","ExomeSeq","ManMut","CNV","RPPA","WB","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.7)
legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","COSMIC","CNV","RPPA","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.6)
legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp","Non-malignant","Unknown"), fill=c("black","red","green","blue","yellow","pink"), title="Subtype", border="white",  bty="n",cex=0.6)
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
legend("topright", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp","Non-malignant","Unknown"), fill=c("black","red","green","blue","yellow","pink"), cex=0.7)
text(x=cutoff_line+1, y=ymax-1, labels=paste("cutoff =", format(mean_cutoff, digits=4), sep=" "), pos=4)
#dev.off()

dev.off()

#Write summary table to file
write.table(predictor_perf, file=results_path, sep="\t", row.names=FALSE)

