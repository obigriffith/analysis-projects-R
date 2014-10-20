#Survival Analysis of cell line predictor for tamoxifen-sensitivity in tamoxifen treated meta-dataset
library(randomForest)
library(survival)

#Set working directory and filenames for Input/output
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/v2/rma_customCDF/")
setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/v2/gcrma_customCDF/")

#Output files
case_pred_outfile="Tamox_CasePredictions.txt"
survival_outfile="Survival.pdf"
survival_outfile2="Survival2.pdf"

#Clinical data
clindatafile="C:/Users/Obi/Documents/My Dropbox/Oana_Obi/Systemically_untreated/Combined_Data/tamoxifen_only-expanded_combined_data.clean.txt"
#clindatafile="C:/Users/Obi/Documents/My Dropbox/Oana_Obi/Systemically_untreated/Combined_Data/tamoxifen_or_endocrine_only-expanded_combined_data.txt"

#Array meta-data (GCRMA normalized, standard CDF - had best estimated performance in cell lines)
#From: /csb/home/oenache/standard_CDF/tamoxifen_only/processing/ALL_gcrma.txt
#datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/ALL_customCDF_tamoxifen-only_rma.txt" #(RMA/customCDF)
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/ALL_gcrma_customCDF.txt" #(GCRMA/customCDF)

#RF Cell line model for tamoxifen sensitivity
#RF_model_file="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/RMA_customCDF/RFmodels/RF_model_Tamoxifen"
RF_model_file="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/drug_response_predictors/U133A/allDrugs/v2/GCRMA_customCDF/RFmodels/RF_model_Tamoxifen"

#File with stringent set of genes used in model building
#stringentgene_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyRMAcustom_genelevel_stringent.csv"
stringentgene_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyGCRMAcustom_genelevel_stringent.csv"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)
cbind(as.vector(clindata[,"GSM"]),header[4:length(header)])

#Get predictor variables
#Selection of stringent set of genes used in model building
stringent_data_import=read.csv(stringentgene_file, row.names=1)
stringent_genes=rownames(stringent_data_import)
stringent_genes_ind=which(rawdata[,3] %in% stringent_genes)
data_stringent=rawdata[stringent_genes_ind,4:length(header)]
rownames(data_stringent)=stringent_genes
predictor_data=t(data_stringent)

#Fix predictor names so that they follow same convention used in rf_model
#colnames(predictor_data)=paste("U133A__",colnames(predictor_data),sep="")

#Load RandomForests classifier from file (object "rf_output" which was saved previously)
load(file=RF_model_file)

#Check that variables are same/overlapping between model and input data
model_variables=names(rf_model$importance[,1])
predictor_variables=colnames(predictor_data)
paste (length(which(model_variables %in% predictor_variables)),"/",length(model_variables),"model variables found in data", sep=" ")

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
clindata_plusRF[tert1,"resistance_group_tert"]=1 #low
clindata_plusRF[tert2,"resistance_group_tert"]=2 #int
clindata_plusRF[tert3,"resistance_group_tert"]=3 #high

#write results to file
write.table(clindata_plusRF,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Create survival plots
#Create new dataframe with just necessary data
tamoxsens_data=clindata_plusRF[,c("general_t_rfs","general_e_rfs","RF_predictions_responses","resistance_group_tert")]

#Censor events beyond 10 yrs
tamoxsens_data[tamoxsens_data[,"general_t_rfs"]>10,"general_e_rfs"]=0
#tamoxsens_data=tamoxsens_data[tamoxsens_data[,"general_t_rfs"]<=10,]

#create a survival object using kidrecurr data
tamoxsens_data.surv = with(tamoxsens_data, Surv(general_t_rfs, general_e_rfs==1))

#Plot KM curve - RF resistance call (binary)
#create kaplan-meier estimate, using gender as strata
krfit.bytamoxsens = survfit(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data)

#Calculate p-value
survdifftest=survdiff(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(survpvalue, digits=4)

#Other tests to try
#coxph(Surv(general_t_rfs, general_e_rfs)~RF_predictions_responses, data=tamoxsens_data, method="breslow")
#survdiff(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data, rho=1)

colors = rainbow(2)

#Manually set Survival p-value to that obtained in SPSS by Breslow test
survpvalue_man=0.063

#plot K-M estimate
pdf(file=survival_outfile)
plot(krfit.bytamoxsens, col = colors, main = "Kaplan-Meier estimate", xlab = "Time (years)", ylab = "Relapse Free Survival")
abline(v = 10, col = "black", lty = 3)
legend(x = "bottomleft", legend = c(levels(tamoxsens_data[,"RF_predictions_responses"]),paste("p=", survpvalue, sep=" ")), col = c(colors,"white"), lty = "solid")
#legend(x = "bottomleft", legend = c(levels(tamoxsens_data[,"RF_predictions_responses"]),paste("p=", survpvalue_man, sep=" ")), col = c(colors,"white"), lty = "solid")
dev.off()


#Plot KM curve - RF resistance group (tertiles)
#create kaplan-meier estimate, using gender as strata
krfit.bytamoxsens = survfit(tamoxsens_data.surv ~ resistance_group_tert, data = tamoxsens_data)

#Calculate p-value
survdifftest=survdiff(tamoxsens_data.surv ~ resistance_group_tert, data = tamoxsens_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(survpvalue, digits=4)
trend_test.by_resistance = coxph(Surv(general_t_rfs, general_e_rfs) ~ resistance_group_tert, data = tamoxsens_data)

colors = rainbow(3)

#plot K-M estimate
pdf(file=survival_outfile2)
plot(krfit.bytamoxsens, col = colors, main = "Kaplan-Meier estimate", xlab = "time (years)", ylab = "Relapse Free Survival")
abline(v = 10, col = "black", lty = 3)
legend(x = "bottomleft", legend = c("low","int","high",paste("p=", survpvalue, sep=" ")), col = c(colors,"white"), lty = "solid")
dev.off()

