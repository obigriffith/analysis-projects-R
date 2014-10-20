library(randomForest)
library(ROCR)
require(Hmisc)

#Set working directory and filenames for Input/output
setwd("/home/obig/Projects/Hodgkins_Lymphoma/")
datafile="expression.100.txt"
performance_outfile="expression.100.RFoutput.txt"
varimp_outfile="expression.100.varImps.txt"
case_pred_outfile="expression.100.CasePredictions.txt"
varimp_pdffile="expression.100.varImps.pdf"
MDS_pdffile="expression.100.MDS.pdf"
case_margins_file="expression.100.Margins.pdf"
ROC_pdffile="expression.100.ROC.pdf"
vote_dist_pdffile="expression.100.vote_dist.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
#To increase speed and reduce memory usage set nrows, colClasses and comment.char as described
#http://www.biostat.jhsph.edu/~rpeng/docs/R-large-tables.html

data=read.table(datafile, header = TRUE, sep="\t", row.names=1, comment.char="")
case_names=colnames(data)
predictor_names=rownames(data)
data=t(data)

#Get potential predictor variables
predictor_data=data[,2:length(data[1,])]

#Get target variable and specify as factor/categorical
target=(data[,1])
target[target==0]="Poor"
target[target==1]="Good"
target=as.factor(target)

#Run RandomForests with imputed marker data as potential predictors and with path_gp1 as target
#rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10000, proximity=TRUE, sampsize=c(17,17))
rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10000, proximity=TRUE, sampsize=c(10,17))

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
write.table(rf_importances[,4],file=varimp_outfile, sep="\t", quote=FALSE, col.names=FALSE)

write("confusion table", file=performance_outfile)
write.table(confusion_out,file=performance_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,err_out,misclass_1,misclass_2), file=performance_outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="Poor"]="P"
target_labels[target_labels=="Good"]="G"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"))
dev.off()

#Determine which cases were misclassified based on OOB testing and by how much
#Margins represent the fraction of votes for the correct class minus the fraction for the incorrect class
#If class is determined by majority rules, then cases with positive margins are correct and vice versa
margins=margin(rf_output, target)
#combine the patient ids with target/known class, predicted class, votes, and margins
case_predictions=cbind(as.vector(target),as.vector(rf_output$predicted),rf_output$votes,as.vector(margins))
colnames(case_predictions)=c("Target","Predicted","Good_votes","Poor_votes","Margins")
misclass_cases=case_predictions[case_predictions[,"Target"]!=case_predictions[,"Predicted"],]

#Write case predictions to file
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Produce graph of margins for all cases
#Note, when plot is called on a 'margin' object, it uses brewer.pal for color selection.
pdf(file=case_margins_file)
plot(margins, ylab="margins (fraction votes true - fraction votes false)", xlab="patient cases")
legend(50,0, legend=c("Good","Poor"), pch=c(20,20), col=c(brewer.pal(3,"Set1")[1], brewer.pal(3,"Set1")[2]))
dev.off()

#Produce back-to-back histogram of vote distributions for Good and Poor
options(digits=2)
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"Poor"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for 100 Hodgkin Lymphomas classified by RF', axes=TRUE, ylab="Fraction votes (Poor)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE)
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE)
dev.off()

#Create ROC curve plot and calculate AUC
#Can use Poor vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling poor vs good
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC)
text(0.6,0.2,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()
