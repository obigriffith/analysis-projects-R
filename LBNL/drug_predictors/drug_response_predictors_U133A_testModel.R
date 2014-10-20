library(randomForest)

#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A_project/Results/")
RF_model_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A_project/Models/RFModel_RMA_standardCDF_AKT12inhibitor"
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered/Neve_AffyRMA_genelevel_maxvar_stringent.csv"

outfile="Cepheid_RFoutput.txt"
case_pred_outfile="Cepheid_CasePredictions.txt"

#Import predictor data - expecting gene names are rows, patient ID as columns
data_import=read.csv(datafile, row.names=1)

#Get one "patient" to simulate single CEL file situation
testdata=as.data.frame(data_import[,1])
rownames(testdata)=rownames(data_import)
colnames(testdata)=colnames(data_import)[1]

#Get predictor variables
predictor_data=t(testdata)

#Load RandomForests classifier from file (object "rf_output" which was saved previously)
load(file=RF_model_file)

#Run test data through forest
RF_predictions_response=predict(rf_model, predictor_data, type="response")
RF_predictions_prob=predict(rf_model, predictor_data, type="prob")
RF_predictions_vote=predict(rf_model, predictor_data, type="vote", norm.votes=FALSE)

#Express results as sensitivity to drug (TRUE/FALSE) and probability of sensitivity to drug
drug_sensitivity=RF_predictions_response=="sensitive"
drug_sensitivity_prob=RF_predictions_prob[1,"sensitive"]
