library(randomForest)
library(ROCR)
require(Hmisc)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimFociNoDiff_vs_DTC")
datafile="ATCprimFociNoDiff_vs_DTC.txt"
outfile="ATCprimFociNoDiff_vs_DTC_RFoutput_all.txt"
outfile="ATCprimFociNoDiff_vs_DTC_RFoutput_transf_markers.txt"
outfile="ATCprimFociNoDiff_vs_DTC_RFoutput_transf_markers_minusVEGF_BCL2.txt"


case_pred_outfile="ATCprimFociNoDiff_vs_DTC_CasePredictions.txt"
varimp_pdffile="ATCprimFociNoDiff_vs_DTC_varImps.pdf"
MDS_pdffile="ATCprimFociNoDiff_vs_DTC_MDS.pdf"
ROC_pdffile="ATCprimFociNoDiff_vs_DTC_ROC.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/RF_analysis")


#Get potential predictor variables (all markers)
predictor_data=data[,10:63]
predictor_names=colnames(predictor_data)

#Get potential predictor variables (just transformation markers)
predictor_data=data[,c("BCL2","ECAD","BetaCatenin","Thyroglobulin","MIB1","P53","VEGF","TOPO2")]
predictor_names=c("Bcl-2","E-CAD","CTNNB1","TG","MIB-1","P53","VEGF","TOPO-II")

#Get potential predictor variables (just transformation markers, minus VEGF and BCL2)
predictor_data=data[,c("ECAD","BetaCatenin","Thyroglobulin","MIB1","P53","TOPO2")]
predictor_names=c("E-CAD","CTNNB1","TG","MIB-1","P53","TOPO-II")

colnames(predictor_data)=predictor_names

#Get target variable and specify as factor/categorical
target=(data[,7])
target[target==0]="ATC"
target[target==1]="DTC"
target=as.factor(target)

#If there are predictor variables that are constant/invariant, consider removing them
#Make sure predictor data is correctly recognized as categorical and determine which variables are invariant (have only one value)
variant_list=vector(mode = "logical", length = 0)
for(i in 1:length(predictor_data)){
 predictor_data[,i]=as.factor(predictor_data[,i])
 if(nlevels(predictor_data[,i])<2){
 variant_list[i]=FALSE
 }
 if(nlevels(predictor_data[,i])>=2){
 variant_list[i]=TRUE
 }
}
predictor_data_variant=predictor_data[,variant_list]

#If there are missing values these will have to be dealt with somehow. Use rfImpute 
#Warning: rfImpute adds the target/response variable as the first column, therefore remove this column with "[,-1]"
predictor_data_imputed=rfImpute(x=predictor_data_variant, y=target, iter=5)[,-1]
#predictor_data_imputed=rfImpute(x=predictor_data, y=target, iter=5)[,-1]

#Run RandomForests with imputed marker data as potential predictors and with path_gp1 as target
rf_output=randomForest(x=predictor_data_imputed, y=target, importance = TRUE, ntree = 10000, proximity=TRUE)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,err_out,misclass_1,misclass_2), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="Benign"]="B"
target_labels[target_labels=="Malignant"]="M"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"))
dev.off()

#Determine which cases were misclassified based on OOB testing and by how much
#Margins represent the fraction of votes for the correct class minus the fraction for the incorrect class
#If class is determined by majority rules, then cases with positive margins are correct and vice versa
margins=margin(rf_output, target)
#combine the patient ids with target/known class, predicted class, votes, and margins
case_predictions=cbind(data[,1:3],target,rf_output$predicted,rf_output$votes,as.vector(margins))
misclass_cases=case_predictions[case_predictions[,"target"]!=case_predictions[,"rf_output$predicted"],]

#Write case predictions to file
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Create ROC curve plot and calculate AUC
#Can use Malignant vote fractions or Malignant-Benign vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling malignant vs benign
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC)
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()