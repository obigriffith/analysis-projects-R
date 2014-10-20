library(randomForest)
library(ROCR)
require(Hmisc)
library(genefilter)
library(mclust)
library(heatmap.plus)
library(fBasics)

#Set working directory and filenames for Input/output
setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/randomtop8/")
setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/randomtop17/")

datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/train_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)
clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.train.2.txt"

testdatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)
testclindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.test.2.txt"

ROC_pdffile="Cepheid_ROC.pdf"

#Set N number of predictors to choose randomly and I number of iterations
#N=8
N=17
I=1000
#I=10

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")
testdata_import=read.table(testdatafile, header = TRUE, na.strings = "NA", sep="\t")
testclin_data_import=read.table(testclindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

testclin_data_order=order(testclin_data_import[,"GSM"])
testclindata=testclin_data_import[testclin_data_order,]
testdata_order=order(colnames(testdata_import)[4:length(colnames(testdata_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
testrawdata=testdata_import[,c(1:3,testdata_order)] #grab first three columns, and then remaining columns in order determined above

header=colnames(rawdata)
testheader=colnames(testrawdata)

#If there are predictor variables that are constant/invariant, consider removing them
#Preliminary gene filtering
X=rawdata[,4:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=rawdata[filt,] 

#Get potential predictor variables
predictor_data=t(filt_Data[,4:length(header)]) #Filtered
#predictor_data=t(rawdata[,4:length(header)]) #Unfiltered

predictor_names=filt_Data[,1] #Filtered, probe ids
#predictor_names=rawdata[,1] #Unfiltered, probe ids
#predictor_names=filt_Data[,3] #Filtered, gene symbols
#predictor_names=rawdata[,3] #Unfiltered, gene symbols
#predictor_names=paste(filt_Data[,3]," (",filt_Data[,1],")", sep="") #Filtered, gene symbol + probe ids
#predictor_names=paste(rawdata[,3]," (",rawdata[,1],")", sep="") #Unfiltered, gene symbol + probe ids
colnames(predictor_data)=predictor_names

#Perform multiple sampling and get AUC then average
sampledata=matrix(data=NA, nrow=I, ncol=3)
for (i in 1:I){
 #Get random set of N predictor variables
 topNopt_probes=sample(predictor_names,size=N,replace=FALSE)
 topNopt_data=rawdata[rawdata[,1]%in%topNopt_probes,]
 predictor_data_topN=t(topNopt_data[,4:length(header)]) #topN optimized list
 predictor_names_topN=topNopt_data[,1] #probe
 colnames(predictor_data_topN)=predictor_names_topN

 testtopNopt_data=testrawdata[testrawdata[,1]%in%topNopt_probes,]
 testpredictor_data_topN=t(testtopNopt_data[,4:length(testheader)]) #topN optimized list
 testpredictor_names_topN=testtopNopt_data[,1] #probe
 colnames(testpredictor_data_topN)=testpredictor_names_topN

 #Filter down to just cases which have 10yr FU (i.e., exclude NAs)
 cases_10yr = !is.na(clindata[,"X10yr_relapse"])
 clindata_10yr=clindata[cases_10yr,]
 predictor_data_10yr=predictor_data_topN[cases_10yr,]

 testcases_10yr = !is.na(testclindata[,"X10yr_relapse"])
 testclindata_10yr=testclindata[testcases_10yr,]
 testpredictor_data_10yr=testpredictor_data_topN[testcases_10yr,]

 #Get target variable and specify as factor/categorical
 target=clindata_10yr[,"X10yr_relapse"] #recurrences after 10yrs not considered events
 #target=clindata_10yr[,"e_rfs"] #All recurrences (even after 10yrs) considered events
 target[target==0]="NoRelapse"
 target[target==1]="Relapse"
 target=as.factor(target)
 
 testtarget=testclindata_10yr[,"X10yr_relapse"] #recurrences after 10yrs not considered events
 #target=clindata_10yr[,"e_rfs"] #All recurrences (even after 10yrs) considered events
 testtarget[testtarget==0]="NoRelapse"
 testtarget[testtarget==1]="Relapse"
 testtarget=as.factor(testtarget)

 #Run RandomForests
 #NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
 #Having an odd number of trees avoids this issue and makes the model fully deterministic
 #rf_model=randomForest(x=predictor_data_10yr, y=target, importance = TRUE, ntree = 10001, proximity=TRUE)

 #Optional: use down-sampling to attempt to compensate for unequal class-sizes
 #tmp = as.vector(table(target))
 #num_classes = length(tmp)
 #min_size = tmp[order(tmp,decreasing=FALSE)[1]]
 #sampsizes = rep(min_size,num_classes)
 #rf_model=randomForest(x=predictor_data_10yr, y=target, importance = TRUE, ntree = 50001, proximity=TRUE, sampsize=sampsizes) 
 #rf_model=randomForest(x=predictor_data_10yr, y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes) 

 #Create model, No downsampling
 rf_model=randomForest(x=predictor_data_10yr, y=target, importance = TRUE, ntree = 10001, proximity=TRUE)

 #Create ROC curve plot and calculate AUC
 #Can use Malignant vote fractions or Malignant-Benign vote fractions as predictive variable
 #The ROC curve will be generated by stepping up through different thresholds for calling malignant vs benign
 predictions=as.vector(rf_model$votes[,2])
 pred=prediction(predictions,target)
 #First calculate the AUC value
 perf_AUC=performance(pred,"auc")
 AUC=perf_AUC@y.values[[1]]
 AUC_out=paste("AUC=",AUC,sep="")
 #Then, plot the actual ROC curve
 perf_ROC=performance(pred,"tpr","fpr")
 
 #Now apply model to independent test data and calculate AUC for that
 #Run test data through forest
 RF_predictions_responses=predict(rf_model, testpredictor_data_10yr, type="response")
 RF_predictions_probs=predict(rf_model, testpredictor_data_10yr, type="prob")
 RF_predictions_votes=predict(rf_model, testpredictor_data_10yr, type="vote")
  
 testrelapse_scores=RF_predictions_probs[,"Relapse"]
 testpred=prediction(testrelapse_scores,testtarget)
 #First calculate the AUC value
 testperf_AUC=performance(testpred,"auc")
 testAUC=testperf_AUC@y.values[[1]]
 testAUC_out=paste("AUC=",testAUC,sep="")
 #Then, plot the actual ROC curve
 testperf_ROC=performance(testpred,"tpr","fpr")

 sampledata[i,1:3]=c(I,AUC,testAUC) 
 print(paste(i,AUC,testAUC,sep=" "))
}
colnames(sampledata)=c("I","AUC", "testAUC")

write.table(sampledata,file="AUC_random_sample_data.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

