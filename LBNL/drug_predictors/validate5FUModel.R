#Survival Analysis of cell line predictor for 5FU-sensitivity in TFAC treated dataset (MDA133 from Liedtke 2010 paper)
library(randomForest)
library(survival)

#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/validation/5FU_in_combination/")

#Output files
case_pred_outfile="5FU_CasePredictions.txt"
survival_outfile="Survival.pdf"

#Clinical data
clindatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/validation/5FU_in_combination/MDA133CompleteInfo20070319.txt"

#Array data (GCRMA normalized, standard CDF - had best estimated performance in cell lines)
#From: /csb/home/oenache/standard_CDF/MDA133/MDA133_ALL_gcrma.txt
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/validation/5FU_in_combination/MDA133_ALL_gcrma.txt" #(GCRMA/standardCDF)

#RF Cell line model for 5FU sensitivity
RF_model_file="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/GCRMA_standardCDF/RFmodels/RF_model_5-FU"

#File with stringent set of genes used in model building
maxvarprobe_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A_project/ExtraFiles/Neve_AffyGCRMA_gene_maxvarprobe_mapping.txt"
stringentgene_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A_project/ExtraFiles/Neve_AffyGCRMA_genelevel_maxvar_stringent.csv"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")


#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"achipcolname"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)
cbind(colnames(rawdata)[4:length(colnames(rawdata))], as.vector(clindata[,"achipcolname"]))

#Get predictor variables
data=rawdata[,4:length(header)]
rownames(data)=rawdata[,1]

#Selection of maxvar probes and then stringent set of genes used in model building
GeneProbeMapping.gcrma=read.table(maxvarprobe_file)
data_maxvarprobe=data[as.vector(GeneProbeMapping.gcrma[-1,2]),]
rownames(data_maxvarprobe)=as.vector(GeneProbeMapping.gcrma[-1,1])
stringent_data_import=read.csv(stringentgene_file, row.names=1)
stringent_genes=rownames(stringent_data_import)
data_stringent=data_maxvarprobe[stringent_genes,]

predictor_data=t(data_stringent)

#Load RandomForests classifier from file (object "rf_model" which was saved previously)
load(file=RF_model_file)

#Run test data through forest
RF_predictions_responses=predict(rf_model, predictor_data, type="response")
RF_predictions_probs=predict(rf_model, predictor_data, type="prob")
RF_predictions_votes=predict(rf_model, predictor_data, type="vote", norm.votes=FALSE)

#Join predictions with clinical data
clindata_plusRF=cbind(clindata,RF_predictions_responses,RF_predictions_probs,RF_predictions_votes)

#Also add grouping of resistance probabilities by tertile
x=clindata_plusRF[,"resistant"]
tertiles=quantile(x, probs=c(1/3,2/3))
tert_cutoff1=tertiles[1]
tert_cutoff2=tertiles[2]
tert1_count=length(which(x <= tert_cutoff1))
tert2_count=length(which(x > tert_cutoff1 & x < tert_cutoff2))
tert3_count=length(which(x >= tert_cutoff2))

#Extract low, int, high samples
tert1=which(x <= tert_cutoff1)
tert2=which(x > tert_cutoff1 & x < tert_cutoff2)
tert3=which(x >= tert_cutoff2)

#Tertiles grouping
clindata_plusRF[tert1,"resistance_group_tert"]="low"
clindata_plusRF[tert2,"resistance_group_tert"]="int"
clindata_plusRF[tert3,"resistance_group_tert"]="high"

#write results to file
write.table(clindata_plusRF,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

