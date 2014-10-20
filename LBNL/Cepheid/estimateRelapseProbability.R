
#Set working directory and filenames for Input/output
setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop17")

datafile="Cepheid_CasePredictions_combined.txt"
RelapseProbabilityPlotfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/Paper/suppl_files/RFRS/RelapseProbabilityPlot.Rdata"
RelapseProbabilityFitfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/Paper/suppl_files/RFRS/RelapseProbabilityFit.Rdata"

resultfile="RelapseProbPredictions.txt"

#Read in pre-determined RFRS probabilities
patient_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

RFRS_values=patient_data[,"Relapse"]

#Load Relapse probability loess fit to allow current predicted relapse probability to be estimated from RFRS values
load(file=RelapseProbabilityFitfile)

#Retrieve prob estimates for all patients
retrieveProb=function(x){
	if (x<0.05){x=0.05} #Can't estimate from loess curve in smallest bin (0 to 0.05), conservatively set to estimate for 0.05
	predict(fit,x)
}
RelapseProbs=sapply(RFRS_values,retrieveProb)
patient_data=cbind(patient_data,RelapseProbs)

#Calculate fraction of patients predicted to have less than X% chance of relapse
perc_patients_5perc_risk=(length(which(RelapseProbs<=5))/length(RelapseProbs))*100
perc_patients_10perc_risk=(length(which(RelapseProbs<=10))/length(RelapseProbs))*100

#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Keep all NoRelapses
NoRelapseCases=which(is.na(patient_data[,"X10yr_relapse"]) | patient_data[,"X10yr_relapse"]==0)
RelapseCases=which(patient_data[,"X10yr_relapse"]==1)

#Downsample Relapse cases so that they represent only 15% of total cases:  x / (429+x) = 0.15 [solving for x, ~76]
#random_RelapseCases=sample(x=RelapseCases, size=76, replace = FALSE, prob = NULL)
#patient_data_down=patient_data[c(NoRelapseCases,random_RelapseCases),]

#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=3)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=76, replace = FALSE, prob = NULL)
 patient_data_down=patient_data[c(NoRelapseCases,random_RelapseCases),]
 RFRS_values=patient_data_down[,"Relapse"]
 RelapseProbs=sapply(RFRS_values,retrieveProb)
 perc_patients_5perc_risk=(length(which(RelapseProbs<=5))/length(RelapseProbs))*100
 perc_patients_10perc_risk=(length(which(RelapseProbs<=10))/length(RelapseProbs))*100
 downsampledata[i,1:3]=c(perc_patients_5perc_risk,perc_patients_10perc_risk,length(RelapseProbs))
}
colnames(downsampledata)=c("5_perc","10_perc","N")
write.table(downsampledata,file=resultfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)




