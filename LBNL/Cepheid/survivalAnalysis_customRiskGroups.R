library(survival)

#Start with case predictions
datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/"

setwd(datadir)
case_pred_outfile="Cepheid_CasePredictions_combined.txt"
KMplotfile="Survival_not_downsampled_RFRS_additional_riskGroups.pdf"

#read RF results and clinical data from file
clindata_plusRF=read.table(case_pred_outfile, header = TRUE, na.strings = "NA", sep="\t")

#Create new risk grouping with additional groups
#Add column for new grouping
quantiles=quantile(clindata_plusRF[,"Relapse"], probs=c(0.2,0.4,0.6,0.8))
clindata_plusRF[,"RF_Group2"]=clindata_plusRF[,"Relapse"]
clindata_plusRF[which(clindata_plusRF[,"Relapse"]<=quantiles[1]),"RF_Group2"]="vlow"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[1] & clindata_plusRF[,"Relapse"]<=quantiles[2]),"RF_Group2"]="low"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[2] & clindata_plusRF[,"Relapse"]<=quantiles[3]),"RF_Group2"]="int"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[3] & clindata_plusRF[,"Relapse"]<=quantiles[4]),"RF_Group2"]="high"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[4]),"RF_Group2"]="vhigh"

#Add column of 10yr censored data
clindata_plusRF[,"e_rfs_10yrcens"]=clindata_plusRF[,"e_rfs"]
clindata_plusRF[which(clindata_plusRF[,"t_rfs"]>10),"e_rfs_10yrcens"]=0

#First, perform survival analysis for entire patient cohort without down-sampling
#Create survival plot and statistics
#Calculate logrank survival statistic between groups
#Create new dataframe with just necessary data
surv_data=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_Group2")]

#create a survival object using data
surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))
#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ RF_Group2, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(as.numeric(survpvalue), digits=3)

#Linear test p-value 
#Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
#See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
#recode  risk groups as 1,2,3
surv_data_lin=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_Group2")]
surv_data_lin[,"RF_Group2"]=as.vector(surv_data_lin[,"RF_Group2"])
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="vlow"),"RF_Group2"]=1
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="low"),"RF_Group2"]=2
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="int"),"RF_Group2"]=3
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="high"),"RF_Group2"]=4
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="vhigh"),"RF_Group2"]=5
surv_data_lin[,"RF_Group2"]=as.numeric(surv_data_lin[,"RF_Group2"])
survpvalue_linear=summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~RF_Group2, data=surv_data_lin))$sctest[3]
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3)

##Plot KM curve
krfit.by_RFgroup = survfit(surv_data.surv ~ RF_Group2, data = surv_data)
pdf(file=KMplotfile)
colors = rainbow(5)
title="Survival by RFRS - Not downsampled"
plot(krfit.by_RFgroup, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
#Set order of categories, categories are by default assigned colors alphabetically by survfit
groups=sort(unique(surv_data[,"RF_Group2"])) #returns unique factor levels sorted alphabetically
names(colors)=groups
groups_custom=c("vlow","low","int","high","vhigh")
colors_custom=colors[groups_custom]
group_sizes_custom=table(surv_data[,"RF_Group2"])[groups_custom]
groups_custom=c("Very Low","Low","Int","High","Very High") #Reset names for consistency with manuscript
legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue_linear, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors_custom,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()

