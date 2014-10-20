#Survival Analysis of cell line predictor for tamoxifen-sensitivity in tamoxifen treated meta-dataset
library(survival)

#Set working directory and filenames for Input/output
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/")
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/endocrine_therapy/")
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/v2/")
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/RMA_stdCDF/top25/")
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/RMA_stdCDF/top228/")
#setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/GCRMA_stdCDF/top27/")
setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/GCRMA_stdCDF/top174/")

#Output files
case_pred_outfile="Tamox_CasePredictions.txt"
survival_outfile="Survival.pdf"
survival_outfile2="Survival2.pdf"
survival_outfile3="Survival_ERBB2neg.pdf"


#Clinical data
#clindatafile="C:/Users/Obi/Documents/My Dropbox/Oana_Obi/Systemically_untreated/Combined_Data/tamoxifen_only-expanded_combined_data.clean.txt"
#clindatafile="C:/Users/Obi/Documents/My Dropbox/Oana_Obi/Systemically_untreated/Combined_Data/tamoxifen_or_endocrine_only-expanded_combined_data.txt"
clindatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/tamoxifen_only-expanded_combined_data.clean.2.txt"

#LSSVM predictions
#LSSVMfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/RMA_stdCDF/top25/LSSVMPredictions_RMAst_TMX_25.csv"
#LSSVMfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/RMA_stdCDF/top228/LSSVMPredictions_RMAst_TMX_228.csv"
#LSSVMfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/GCRMA_stdCDF/top27/LSSVMPredictions_GCRMAst_TMX_27.csv"
LSSVMfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TumorData/tamoxifen/LSSVM/GCRMA_stdCDF/top174/LSSVMPredictions_GCRMAst_TMX_174.csv"

#Read in data (expecting a tab-delimited file with header line and rownames)
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")
LSSVM_import=read.csv(LSSVMfile)

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
LSSVM_order=order(LSSVM_import[,"TMX.Sample"])
LSSVM_data=LSSVM_import[LSSVM_order,]
cbind(as.vector(LSSVM_data[,"TMX.Sample"]),as.vector(clindata[,"GSM"]))

#Join predictions with clinical data
clindata_plusLSSVM=cbind(clindata,LSSVM_data)


#Also add grouping of resistance probabilities by tertile
x=clindata_plusLSSVM[,"Predicted.probability.of.response"] #Note probability of reponse is opposite of probability of resistance
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
clindata_plusLSSVM[tert1,"resistance_group_tert"]=1 #high
clindata_plusLSSVM[tert2,"resistance_group_tert"]=2 #int
clindata_plusLSSVM[tert3,"resistance_group_tert"]=3 #low

#Recode LSSVM prediction categories
clindata_plusLSSVM[,"LSSVM_predictions_responses"]=clindata_plusLSSVM[,"Sens.1..Res.0."]
clindata_plusLSSVM[which(clindata_plusLSSVM[,"LSSVM_predictions_responses"]==1),"LSSVM_predictions_responses"]="sensitive"
clindata_plusLSSVM[which(clindata_plusLSSVM[,"LSSVM_predictions_responses"]==0),"LSSVM_predictions_responses"]="resistant"

#write results to file
write.table(clindata_plusLSSVM,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#Create survival plots
#Create new dataframe with just necessary data
tamoxsens_data=clindata_plusLSSVM[,c("general_t_rfs","general_e_rfs","LSSVM_predictions_responses","resistance_group_tert")]
tamoxsens_data_ERBB2neg=clindata_plusLSSVM[which(clindata_plusLSSVM[,"ERBB2amp_array"]==0),c("general_t_rfs","general_e_rfs","LSSVM_predictions_responses","resistance_group_tert")]

#Censor events beyond 10 yrs
tamoxsens_data[tamoxsens_data[,"general_t_rfs"]>10,"general_e_rfs"]=0
tamoxsens_data_ERBB2neg[tamoxsens_data_ERBB2neg[,"general_t_rfs"]>10,"general_e_rfs"]=0
#tamoxsens_data=tamoxsens_data[tamoxsens_data[,"general_t_rfs"]<=10,]

#create a survival object using data
tamoxsens_data.surv = with(tamoxsens_data, Surv(general_t_rfs, general_e_rfs==1))
tamoxsens_data_ERBB2neg.surv = with(tamoxsens_data_ERBB2neg, Surv(general_t_rfs, general_e_rfs==1))

#Plot KM curve - RF resistance call (binary)
#create kaplan-meier estimate, using gender as strata
krfit.bytamoxsens = survfit(tamoxsens_data.surv ~ LSSVM_predictions_responses, data = tamoxsens_data)
krfit.bytamoxsens_ERBB2neg = survfit(tamoxsens_data_ERBB2neg.surv ~ LSSVM_predictions_responses, data = tamoxsens_data_ERBB2neg)

#Calculate p-value
survdifftest=survdiff(tamoxsens_data.surv ~ LSSVM_predictions_responses, data = tamoxsens_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

survdifftest_ERBB2neg=survdiff(tamoxsens_data_ERBB2neg.surv ~ LSSVM_predictions_responses, data = tamoxsens_data_ERBB2neg)
survpvalue_ERBB2neg = 1 - pchisq(survdifftest_ERBB2neg$chisq, length(survdifftest_ERBB2neg$n) - 1)
survpvalue_ERBB2neg = round(as.numeric(survpvalue_ERBB2neg), digits=3)

#Other tests to try
#coxph(Surv(general_t_rfs, general_e_rfs)~RF_predictions_responses, data=tamoxsens_data, method="breslow")
#survdiff(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data, rho=1)

colors = rainbow(2)

#Manually set Survival p-value to that obtained in SPSS by Breslow test
#survpvalue_man=0.063

#plot K-M estimate
pdf(file=survival_outfile)
plot(krfit.bytamoxsens, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
legend(x = "bottomleft", legend = c(unique(tamoxsens_data[,"LSSVM_predictions_responses"]),paste("p=", survpvalue, sep=" ")), col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)
#legend(x = "bottomleft", legend = c(levels(tamoxsens_data[,"LSSVM_predictions_responses"]),paste("p=", survpvalue_man, sep=" ")), col = c(colors,"white"), lty = "solid")
dev.off()

#plot K-M estimate
pdf(file=survival_outfile3)
plot(krfit.bytamoxsens_ERBB2neg, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
legend(x = "bottomleft", legend = c(unique(tamoxsens_data_ERBB2neg[,"LSSVM_predictions_responses"]),paste("p=", survpvalue_ERBB2neg, sep=" ")), col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)
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
legend(x = "bottomleft", legend = c("high","int","low",paste("p=", survpvalue, sep=" ")), col = c(colors,"white"), lty = "solid")
dev.off()





