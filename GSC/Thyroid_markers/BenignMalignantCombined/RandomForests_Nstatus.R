library(randomForest)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"
outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_Nstatus_RFoutput.txt"
varimp_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_Nstatus_varImps.pdf"
MDS_pdffile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_Nstatus_MDS.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
data=data[data[,13]<4,] #Get rid of M samples
data=data[data[,13]>0,] #Also exclude benign

#Get potential predictor variables
predictor_data=data[,36:92]

#Get target variable and specify as factor/categorical
target=(data[,"N"])
target[target==0]="N0"
target[target==1]="N1"
target=as.factor(target)

#Only use data where target is not NA
predictor_data=predictor_data[complete.cases(target),]
target=target[complete.cases(target)]

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

#If there are missing values these will have to be dealt with somehow.  
#Warning: rfImpute adds the target/response variable as the first column, therefore remove this column with "[,-1]"
#predictor_data_imputed=rfImpute(x=predictor_data_variant, y=target, iter=5)[,-1]
predictor_data_imputed=rfImpute(x=predictor_data, y=target, iter=5)[,-1]

#If there are predictor variables that are constanct/invariant, consider removing them

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

#Produce graph of variable importances
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="N0"]="N0"
target_labels[target_labels=="N1"]="N1"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"))
dev.off()
