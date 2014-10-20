#This script will summarize the concordance between drug response predictions versus measured sensitivity for ~10 independant lines

#Load drug response metrics
drugdata_file="/Users/ogriffit/Dropbox/drug_predictors/validation/newdrugresponsedata.txt"
drugdata=read.table(file=drugdata_file, header=TRUE, sep="\t", row.names=1)

#Load GI50 values to assign classes
meanGI50data_file="/Users/ogriffit/Dropbox/drug_predictors/validation/GI50meanThresholds_v2.csv"
meanGI50data=read.csv(file=meanGI50data_file, row.names=1)

#Set working dir for results
setwd("/Users/ogriffit/Dropbox/drug_predictors/validation/RtoolboxResults/Breast/Results/")

#Load drug response predictions from Gray lab Rtoolbox
Rtoolbox_results_file="validation_results_ExpCNV.txt"
predictions=read.table(file=Rtoolbox_results_file, header=TRUE, sep="\t")

#Drugs to analyze for concordance
#drugs=c("API-2(Tricir)","BIBW2992","GSK2126458A","GSK2141795c","GSK_Tykerb","Iressa","Rapamycin")
drugs=c("BIBW2992","GSK2126458A","GSK2141795c","GSK_Tykerb","Iressa","Rapamycin") #Exclude Tricir since its best model has AUC < 0.6
drugdata=drugdata[,drugs]

#Cell lines to analyze
cell_lines=c("21MT1","21MT2","21NT","21PT","AU565","EFM192A","EFM192B","EFM192C","HCC202","HCC3153","JIMT1")

#Set up dataframe to store results
#Keep existing column names and datatypes 
results=predictions[1:length(which(!is.na(drugdata))),]
#Add extra columns for additional information
results=cbind(results,matrix(NA, nrow = length(which(!is.na(drugdata))), ncol=2))
results[1:length(rownames(results)),1:length(colnames(results))]=NA
colnames(results)=c(colnames(predictions),"GI50","meanGI50")

k=0
for (i in 1:length(drugs)){
 drug=drugs[i]
 for (j in 1:length(cell_lines)){
  cell_line=cell_lines[j]
  if (!is.na(drugdata[cell_line,drug])){
   k=k+1
   GI50=drugdata[cell_line,drug]
   meanGI50=meanGI50data[drug,1]
   prediction=predictions[which(predictions[,"Sample"]==cell_line & predictions[,"Drug.Compound"]==drug),]
   results[k,]=c(prediction,GI50,meanGI50)
  }	
 }
}

#Determine measured sensitivity by comparing GI50 values to mean GI50 cutoffs
resistants=which(results[,"GI50"]<=results[,"meanGI50"])
sensitives=which(results[,"GI50"]>results[,"meanGI50"])
results[,"MeasuredSensitivity"]="NA"
results[resistants,"MeasuredSensitivity"]="resistant"
results[sensitives,"MeasuredSensitivity"]="sensitive"

#For ease of interpretation, convert predicted sensitivity +/- to words
results[,"Sensitivity"]=as.vector(results[,"Sensitivity"])

results[which(results[,"Sensitivity"]=="+"),"Sensitivity"]="sensitive"
results[which(results[,"Sensitivity"]=="-"),"Sensitivity"]="resistant"

table(results[,c("Sensitivity","MeasuredSensitivity")])
write.table(results, file="validation_results.txt",row.names=FALSE, col.names=TRUE, sep="\t")


