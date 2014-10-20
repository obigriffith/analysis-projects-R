library(survival)

#Start with case predictions
#datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/test_survival/allgene/" #1
#datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/test_survival/17gene_optimized/" #2
#datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/test_survival/8gene_optimized/" #3
datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/" #4

#Set type for either train or test. Different downsampling command below for each.
#type="test" #1,2,3
type="train" #4

setwd(datadir)
#case_pred_outfile="Cepheid_CasePredictions.txt" #1,2,3
case_pred_outfile="Cepheid_CasePredictions_combined.txt" #4
#case_pred_1000downsamp_10yrRFS_KMpvalue_outfile="Cepheid_CasePredictions_1000downsamp_10yrRFS_KMpvalue.txt" #1,2,3
case_pred_1000downsamp_10yrRFS_KMpvalue_outfile="Cepheid_CasePredictions_combined_1000downsamp_10yrRFS_KMpvalue.txt" #4
KMplotfile="Survival_not_downsampled_RFRS.pdf"

#read RF results and clinical data from file
clindata_plusRF=read.table(case_pred_outfile, header = TRUE, na.strings = "NA", sep="\t")

#Fix column names from older files
colnames(clindata_plusRF)[which(colnames(clindata_plusRF)=="RF_Group1")]="RF_risk_group"

#Add column of 10yr censored data
clindata_plusRF[,"e_rfs_10yrcens"]=clindata_plusRF[,"e_rfs"]
clindata_plusRF[which(clindata_plusRF[,"t_rfs"]>10),"e_rfs_10yrcens"]=0

#First, perform survival analysis for entire patient cohort without down-sampling
#Create survival plot and statistics
#Calculate logrank survival statistic between groups
#Create new dataframe with just necessary data
surv_data=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_risk_group")]

#create a survival object using data
surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))
#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ RF_risk_group, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(as.numeric(survpvalue), digits=3)

#Linear test p-value 
#Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
#See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
#recode  risk groups as 1,2,3
surv_data_lin=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_risk_group")]
surv_data_lin[,"RF_risk_group"]=as.vector(surv_data_lin[,"RF_risk_group"])
surv_data_lin[which(surv_data_lin[,"RF_risk_group"]=="low"),"RF_risk_group"]=1
surv_data_lin[which(surv_data_lin[,"RF_risk_group"]=="int"),"RF_risk_group"]=2
surv_data_lin[which(surv_data_lin[,"RF_risk_group"]=="high"),"RF_risk_group"]=3
surv_data_lin[,"RF_risk_group"]=as.numeric(surv_data_lin[,"RF_risk_group"])
survpvalue_linear=summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~RF_risk_group, data=surv_data_lin))$sctest[3]
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3)

##Plot KM curve
krfit.by_RFgroup = survfit(surv_data.surv ~ RF_risk_group, data = surv_data)
pdf(file=KMplotfile)
colors = rainbow(3)
title="Survival by RFRS - Not downsampled"
plot(krfit.by_RFgroup, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
#Set order of categories, categories are by default assigned colors alphabetically by survfit
groups=levels(surv_data[,"RF_risk_group"]) #returns unique factor levels sorted alphabetically
names(colors)=groups
groups_custom=c("low","int","high")
colors_custom=colors[groups_custom]
group_sizes_custom=table(surv_data[,"RF_risk_group"])[groups_custom]
groups_custom=c("Low","Int","High") #Reset names for consistency with manuscript
legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue_linear, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors_custom,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()







#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Keep all NoRelapses
NoRelapseCases=which(is.na(clindata_plusRF[,"X10yr_relapse"]) | clindata_plusRF[,"X10yr_relapse"]==0)
RelapseCases=which(clindata_plusRF[,"X10yr_relapse"]==1)

#Downsample Relapse cases so that they represent only 15% of total cases (all non-10yr_relapse):  x / (216+x) = 0.15 [solving for x, ~38]
#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
#I=10 #for testing
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=8)
for (i in 1:I){
 if (type=="test"){
   random_RelapseCases=sample(x=RelapseCases, size=38, replace = FALSE, prob = NULL)
 }
 if (type=="train"){
   random_RelapseCases=sample(x=RelapseCases, size=76, replace = FALSE, prob = NULL)
 }
 case_predictions_all_combined_down=clindata_plusRF[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="low","X10yr_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="int","X10yr_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group"]=="high","X10yr_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100

 #Calculate logrank survival statistic between groups
 #Create new dataframe with just necessary data
 surv_data=case_predictions_all_combined_down[,c("t_rfs","e_rfs_10yrcens","RF_risk_group")]
 #create a survival object using data
 surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))
 #Calculate p-value
 survdifftest=survdiff(surv_data.surv ~ RF_risk_group, data = surv_data)
 survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)

 #Linear test p-value 
 #Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
 #See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
 #recode  risk groups as 1,2,3
 surv_data[,"RF_risk_group"]=as.vector(surv_data[,"RF_risk_group"])
 surv_data[which(surv_data[,"RF_risk_group"]=="low"),"RF_risk_group"]=1
 surv_data[which(surv_data[,"RF_risk_group"]=="int"),"RF_risk_group"]=2
 surv_data[which(surv_data[,"RF_risk_group"]=="high"),"RF_risk_group"]=3
 surv_data[,"RF_risk_group"]=as.numeric(surv_data[,"RF_risk_group"])
 survpvalue_linear=summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~RF_risk_group, data=surv_data))$sctest[3]

 #survpvalue = round(as.numeric(survpvalue), digits=3)
 ##Plot KM curve
 #krfit.by_RFgroup = survfit(surv_data.surv ~ RF_risk_group, data = surv_data)
 #colors = rainbow(3)
 #title="Survival by RFRS"
 #plot(krfit.by_RFgroup, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
 #abline(v = 10, col = "black", lty = 3)
 #groups=unique(surv_data[,"RF_risk_group"])
 #group_sizes=table(surv_data[,"RF_risk_group"])[groups]
 #legend_text=c(paste(groups, " ", "(", group_sizes, ")", sep=""),paste("p =", survpvalue, sep=" "))
 #legend(x = "bottomleft", legend = legend_text, col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)

 downsampledata[i,1:8]=c(perc_10yr_relapse_low,length(low_10yr_relapses),perc_10yr_relapse_int,length(int_10yr_relapses),perc_10yr_relapse_high,length(high_10yr_relapses),survpvalue,survpvalue_linear)
}
colnames(downsampledata)=c("low_perc","low_N","int_perc","int_N","high_perc","high_N","p-value","p-value(linear_trend)")

#Print means to screen
apply(downsampledata,2,mean)

#Print results to file
write.table(downsampledata,file=case_pred_1000downsamp_10yrRFS_KMpvalue_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Other survival tests to try
#coxph(Surv(general_t_rfs, general_e_rfs)~RF_predictions_responses, data=tamoxsens_data, method="breslow") #linear trend test?
#survdiff(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data, rho=1)

