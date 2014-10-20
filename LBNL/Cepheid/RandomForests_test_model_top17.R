library(randomForest)
library(ROCR)
require(Hmisc)

#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/test_survival/17gene_optimized/")

RF_model_file="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop17/RF_model_17gene_optimized"

datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/processing/processed_final2/test_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)
clindatafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.test.2.txt"

outfile="Cepheid_RFoutput.txt"
case_pred_outfile="Cepheid_CasePredictions.txt"
case_margins_file="Cepheid_Margins.pdf"
ROC_pdffile="Cepheid_ROC.pdf"
vote_dist_pdffile="Cepheid_vote_dist.pdf"
combined_case_pred_downsamp_outfile="Cepheid_CasePredictions_combined_downsamp.txt"
combined_case_pred_1000downsamp_10yrRFS_outfile="Cepheid_CasePredictions_combined_1000downsamp_10yrRFS.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)

#Get predictor variables
top17opt_probes=c("204767_s_at","10682_at","201291_s_at","9133_at","1164_at","208079_s_at","23224_at","55435_at","23220_at","201461_s_at","202709_at","57122_at","23405_at","201483_s_at","29127_at","204416_x_at","10628_at")
top17opt_data=rawdata[rawdata[,1]%in%top17opt_probes,]
predictor_data=t(top17opt_data[,4:length(header)]) #Top20 optimized list
predictor_names=top17opt_data[,3] #gene symbol
colnames(predictor_data)=predictor_names

#Load RandomForests classifier from file (object "rf_model" which was saved previously)
load(file=RF_model_file)

#Run test data through forest
RF_predictions_responses=predict(rf_model, predictor_data, type="response")
RF_predictions_probs=predict(rf_model, predictor_data, type="prob")
RF_predictions_votes=predict(rf_model, predictor_data, type="vote")

#Determine RF risk score according to previously determined thresholds
#low < 0.333
#0.333 >= int < 0.606 
#high >= 0.606

RF_risk_group=RF_predictions_votes[,"Relapse"]
RF_risk_group[RF_predictions_votes[,"Relapse"]<0.333]="low"
RF_risk_group[RF_predictions_votes[,"Relapse"]>=0.333 & RF_predictions_votes[,"Relapse"]<0.606]="int"
RF_risk_group[RF_predictions_votes[,"Relapse"]>=0.606]="high"


#Join predictions with clinical data
clindata_plusRF=cbind(clindata,RF_predictions_responses,RF_predictions_votes,RF_risk_group)

#write results to file
write.table(clindata_plusRF,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Keep all NoRelapses
NoRelapseCases=which(is.na(clindata_plusRF[,"X10yr_relapse"]) | clindata_plusRF[,"X10yr_relapse"]==0)
RelapseCases=which(clindata_plusRF[,"X10yr_relapse"]==1)




#Downsample Relapse cases so that they represent only 15% of total cases (all non-10yr_relapse):  x / (216+x) = 0.15 [solving for x, ~38]
#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=6)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=38, replace = FALSE, prob = NULL)
 case_predictions_all_combined_down=clindata_plusRF[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="low","X10yr_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="int","X10yr_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="high","X10yr_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
 downsampledata[i,1:6]=c(perc_10yr_relapse_low,length(low_10yr_relapses),perc_10yr_relapse_int,length(int_10yr_relapses),perc_10yr_relapse_high,length(high_10yr_relapses))
}
colnames(downsampledata)=c("low_perc","low_N","int_perc","int_N","high_perc","high_N")
#Print means to screen
mean(downsampledata[,"low_perc"]);mean(downsampledata[,"int_perc"]);mean(downsampledata[,"high_perc"])
write.table(downsampledata,file=combined_case_pred_1000downsamp_10yrRFS_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Create representative result for plotting survival figure: run the block below until %'s similar to means above
random_RelapseCases=sample(x=RelapseCases, size=38, replace = FALSE, prob = NULL)
case_predictions_all_combined_down=clindata_plusRF[c(NoRelapseCases,random_RelapseCases),]
low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="low","X10yr_relapse"]
int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="int","X10yr_relapse"]
high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="high","X10yr_relapse"]
perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
perc_10yr_relapse_low; perc_10yr_relapse_int; perc_10yr_relapse_high #Check values against means for above
write.table(case_predictions_all_combined_down,file=combined_case_pred_downsamp_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Determine performance statistics
#Use only patients with 10 yr FU

RF_predictions_vs_relapse_10yrFU=clindata_plusRF[!is.na(clindata_plusRF[,"X10yr_relapse"]),c("X10yr_relapse","RF_predictions_responses","Relapse","RF_risk_group")]
confusion=table(RF_predictions_vs_relapse_10yrFU[,1:2])
rownames(confusion)=c("NoRelapse","Relapse")
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=((confusion[1,2]+confusion[2,1])/sum(confusion))*100
overall_accuracy=((confusion[1,1]+confusion[2,2])/sum(confusion))*100
class1_error=confusion[1,2]/(confusion[1,1]+confusion[1,2])
class2_error=confusion[2,1]/(confusion[2,2]+confusion[2,1])

#Get performance stats for individual risk groups
#low risk
RF_predictions_vs_relapse_10yrFU_low=RF_predictions_vs_relapse_10yrFU[RF_predictions_vs_relapse_10yrFU[,"RF_risk_group"]=="low",]
confusion_low=table(RF_predictions_vs_relapse_10yrFU_low[,1:2])
rownames(confusion_low)=c("NoRelapse","Relapse")
overall_error_low=((confusion_low[1,2]+confusion_low[2,1])/sum(confusion_low))*100
overall_accuracy_low=((confusion_low[1,1]+confusion_low[2,2])/sum(confusion_low))*100

#high risk
RF_predictions_vs_relapse_10yrFU_high=RF_predictions_vs_relapse_10yrFU[RF_predictions_vs_relapse_10yrFU[,"RF_risk_group"]=="high",]
confusion_high=table(RF_predictions_vs_relapse_10yrFU_high[,1:2])
rownames(confusion_high)=c("NoRelapse","Relapse")
overall_error_high=((confusion_high[1,2]+confusion_high[2,1])/sum(confusion_high))*100
overall_accuracy_high=((confusion_high[1,1]+confusion_high[2,2])/sum(confusion_high))*100

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
class1_error_out=paste(colnames(confusion)[1]," error rate= ",class1_error, sep="")
class2_error_out=paste(colnames(confusion)[2]," error rate= ",class2_error, sep="")
misclass_1=paste(confusion[1,2], colnames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], colnames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")
err_out_low=paste("overall error rate (low risk group)=",overall_error_low,sep="")
acc_out_low=paste("overall accuracy (low risk group)=",overall_accuracy_low,sep="")
err_out_high=paste("overall error rate (high risk group)=",overall_error_high,sep="")
acc_out_high=paste("overall accuracy (high risk group)=",overall_accuracy_high,sep="")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Produce back-to-back histogram of vote distributions for Relapse and NoRelapse
target=RF_predictions_vs_relapse_10yrFU[,"X10yr_relapse"]
target[target==1]="Relapse"
target[target==0]="NoRelapse"
relapse_scores=RF_predictions_vs_relapse_10yrFU[,"Relapse"]
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(relapse_scores, target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes (Relapse)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()

#Create ROC curve plot and calculate AUC
#Can use Relapse vote fractions fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling Response vs NoResponse

pred=prediction(relapse_scores,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
AUC_out=paste("AUC=",AUC,sep="")
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()

#Print results to file
write("confusion table", file=outfile)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error_out,class2_error_out,misclass_1,misclass_2,AUC_out,err_out_low,acc_out_low,err_out_high,acc_out_high), file=outfile, append=TRUE)

