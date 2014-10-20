library(randomForest)
library(mclust)
require(Hmisc)
library(ROCR)

#drug="private_company_drug_8271" #pan PI3K
drug="private_company_drug_2643" #Lapatinib

outdir=paste("C:/Users/Obi/Documents/Projects/BCCL/drug_response_predictors/",drug,"/",sep="")
setwd(outdir)
outfile=paste("RFoutput.txt_",drug,".txt",sep="")
case_pred_outfile="CasePredictions.txt"
varimp_pdffile="varImps.pdf"
MDS_pdffile="MDS.pdf"
case_margins_file="Margins.pdf"
ROC_pdffile="ROC.pdf"
vote_dist_pdffile="vote_dist.pdf"
vote_density_pdffile="density_RF.pdf"
margins_mixed_model_clustering_pdffile="Margins_MM.pdf"


#Parameters
pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 2 #Minimum coefficient of variation (0.7 recommended?)
cov_max = 10 #Maximum cov

#Libraries that have drug response, RNAseq and mutation data
complete_libs=c("BT20","BT474","BT549","CAMA1","HCC1143","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC202","HCC3153","HCC38","HCC70","HS578T","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SUM1315MO2","SUM149PT","SUM159PT","T47D","UACC812","UACC893","ZR751","ZR7530")

#Libraries that have drug response and RNAseq data
RNAseq_libs=c("X184A1","X184B5","X600MPE","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SKBR3","SUM1315","SUM149PT","SUM159PT","SUM225CWN","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")


#RNAseq data (Choose one feature type)
#Gene
#datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_GeneExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_GeneExpression_v53.txt"

#Transcript
#datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_TranscriptExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_TranscriptExpression_v53.txt"

#Junction
datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_KnownJunctionExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_KnownJunctionExpression_v53.txt"

raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:62]
raw_exp_status=raw_exp_status_import[,5:62]

#Fix misnamed library
colnames(raw_data)[which(colnames(raw_data)=="MDAMB13v1")]="MDAMB134VI"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="MDAMB13v1")]="MDAMB134VI"

#Retrieve data for only libraries with both RNAseq and drug response data
data=raw_data[,RNAseq_libs]
feat_data=raw_feat_data
exp_status=raw_exp_status[,RNAseq_libs]
libs=colnames(data)

#Define a percent expressed function and filter out features with less than minimum
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}
pe_data=apply(exp_status, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data_filt=data[passed_pe,]
feat_data_filt=feat_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
  return(cov)
}

cov_data=apply(data_filt, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data_filt=data_filt[passed_cov,]
feat_data_filt=feat_data_filt[passed_cov,]


#rownames(data_filt)=feat_data_filt[,"FID"]
#rownames(data_filt)=feat_data_filt[,"Seq_Name"]
rownames(data_filt)=paste(feat_data_filt[,"EnsEMBL_Gene_ID"],feat_data_filt[,"Seq_Name"], sep="_")

#drug response data
drugdatafile="C:/Users/Obi/Documents/Projects/BCCL/drug_data/LBNL_Gray_BCCL_GI50_10Feb_final_encoded.txt"
raw_drugdata_import=read.table(drugdatafile, header = TRUE, na.strings = "NA", sep="\t")

drug_data=raw_drugdata_import[,2:75]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Retrieve data for only libraries with both RNAseq and drug response data
drug_data_filt=drug_data[RNAseq_libs,]

#Transform data to -log values
drug_data_filt_trans=-log(drug_data_filt)
drugs=colnames(drug_data_filt_trans)


#For drug of interest, find expression features associated with response and build predictor for drug response
#Divide cell lines into responders and non-responders for drug

x=as.numeric(drug_data_filt_trans[,drug])
x_NA=which(is.na(x))
x_ordered=order(x,na.last=NA)
cutoff=round(length(x_ordered)/2)
resistants=x_ordered[1:cutoff]
sensitives=x_ordered[(cutoff+1):length(x_ordered)]
response_class=vector(length=length(x))
response_class[x_NA]=NA
response_class[sensitives]="sensitive"
response_class[resistants]="resistant"

#Exclude libs where response_class=NA
nonNA=which(!is.na(response_class))
data_filt_nonNA=data_filt[,nonNA]
response_class_nonNA=response_class[nonNA]

#Make sure there are at least 10 libs classified as sensitive and resistant
length(which(response_class_nonNA=="sensitive"))
length(which(response_class_nonNA=="resistant"))

target=as.factor(response_class_nonNA)

#Create RF predictor for response_class using RNAseq data(data_filt)
rf_output=randomForest(x=t(data_filt_nonNA), y=target, importance = TRUE, ntree = 10001, proximity=TRUE)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=1-overall_error
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="sensitive"]="S"
target_labels[target_labels=="resistant"]="R"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
dev.off()


#Produce back-to-back histogram of vote distributions for Sensitive and Resistant
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"sensitive"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes (sensitive)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()


#Create ROC curve plot and calculate AUC
#Can use Sensitive vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling sensitive vs resistant
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
