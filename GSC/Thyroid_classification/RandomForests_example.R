library(randomForest)
library(ROCR)
require(Hmisc)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis/data_files/")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"
outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput.txt"
case_pred_outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_CasePredictions.txt"
varimp_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_varImps.pdf"
MDS_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_MDS.pdf"
case_margins_file="BenignMalignant_58markers_23JAN08_ungrouped_NoM_Margins.pdf"
ROC_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_ROC.pdf"
vote_dist_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_vote_dist.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
data=data[data[,13]<4,] #Keep only rows where pathology is 0, 1, 2 or 3

#Get potential predictor variables
predictor_data=data[,36:92]

predictor_names=c("AAT","AMF-R","AR","Aurora-A","Aurora-C","Bcl-2","CTNNB1","CDX2","CK19","c-kit","COX2","CR3","Cyclin-D1","Cyclin-E","Caveolin","CAV-1","Clusterin","E-CAD","ER","Galectin-3","HBME-1","EGFR","HER2","HER3","HER4","HSP-27","IGFBP2","IGFBP5","INH","KI67","MDM2","MLH1","MRAS","O13","P16","P21","P27","P53","P57","P63","P75-NTR","PGI","PMS2","PR","PSA","RET","S100","Syntrophin","TDT","TG","TOPO-II","TS106","TSH","TTF-1","VEGF","WT1","P504S")
predictor_names=c("Galectin-3","Aurora-A","P16","AR","HBME-1","Bcl-2","Cyclin-D1","CAV-1","Cyclin-E","E-CAD","CR3","Clusterin","IGFBP5","P21","IGFBP2","CTNNB1","HER4","TG","KI67","Caveolin","Aurora-C","S100","MRAS","c-kit","HER3","RET","AMF-R","MLH1","AAT","TTF-1","PGI","HSP-27","Syntrophin")

colnames(predictor_data)=predictor_names

#Get target variable (path_gp1) and specify as factor/categorical
target=(data[,14])
target[target==0]="Benign"
target[target==1]="Malignant"
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

#Produce graph of margins for all cases
#Note, when plot is called on a 'margin' object, it uses brewer.pal for color selection.
pdf(file=case_margins_file)
plot(margins, ylab="margins (fraction votes true - fraction votes false)", xlab="patient cases")
legend(100,0, legend=c("benign","malignant"), pch=c(20,20), col=c(brewer.pal(3,"Set1")[1], brewer.pal(3,"Set1")[2]))
dev.off()

#Produce back-to-back histogram of vote distributions for Benign and Malignant
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"Malignant"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for 199 thryoid lesions classified by RF', axes=TRUE, ylab="Fraction votes (Malignant)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()

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