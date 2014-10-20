library(randomForest)
library(preprocessCore)
library(ROCR)
library(survival)
library(mclust)
library(illuminaHumanv3.db)
library(plotrix)

#Load training data - for normalization purposes
datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
header=colnames(data_import)
training_data=data_import[,4:length(header)]
training_features=data_import[,1:3]

#Load METABRIC test data
METABRIC_clinfeatfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/independent_research/Complete_METABRIC_Clinical_Features_Data.rbin"
load(METABRIC_clinfeatfile)
summary(Complete_METABRIC_Clinical_Features_Data)
METABRIC_DSSfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/independent_research/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin"
load(METABRIC_DSSfile)
summary(Complete_METABRIC_Clinical_Survival_Data__DSS)

METABRIC_OSfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/independent_research/Complete_METABRIC_Clinical_Survival_Data_OS.rbin"
load(METABRIC_OSfile)
summary(Complete_METABRIC_Clinical_Survival_Data_OS)

METABRIC_expfile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/independent_research/Complete_METABRIC_Expression_Data.rbin"
load(METABRIC_expfile)
summary(Complete_METABRIC_Expression_Data)

expression_data=exprs(Complete_METABRIC_Expression_Data)
clinical_data_raw=Complete_METABRIC_Clinical_Features_Data
probe_ids=rownames(expression_data)
probe_mapping=unlist(mget(probe_ids, illuminaHumanv3SYMBOL, ifnotfound=NA))

#Load RandomForests classifier from file (object "rf_model" which was saved previously)
#RF_model_file="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop17/RF_model_17gene_optimized"
RF_model_file="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop8/RF_model_8gene_optimized"
#RF_model_file="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees_repeat/RF_model"
rf_model_object=load(file=RF_model_file)
rf_model=get(rf_model_object)

#Set working dir for results
#result_dir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/results/17gene"
#result_dir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/results/8gene"
#result_dir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/results/17gene_optimized"
result_dir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/METABRIC/results/8gene_optimized"

setwd(result_dir)
case_pred_outfile="METABRIC_CasePredictions.txt"
clindata_file="METABRIC_complete_clinical_data.txt"
ROC_pdffile="METABRIC_ROC.pdf"
KMplotfile="METABRIC_survival.pdf"
case_pred_1000downsamp_KMpvalue_outfile="METABRIC_Cepheid_CasePredictions_1000downsamp_KMpvalue.txt"

#Print out starting point for clinical data to file
write.table(clinical_data_raw,file=clindata_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Need to strip down to just LN-, ER+, HER2- (by probe or amplicon) who did not receive chemo
LNneg_ERpos_HER2neg=which(clinical_data_raw[,"lymph_nodes_positive"]==0 & clinical_data_raw[,"ER_IHC_status"]=="pos" & clinical_data_raw[,"ER.Expr"]=="+" & clinical_data_raw[,"HER2_SNP6_state"]!="GAIN" & clinical_data_raw[,"Her2.Expr"]=="-" & (clinical_data_raw[,"Treatment"]=="HT" | clinical_data_raw[,"Treatment"]=="NONE"))
clinical_data=clinical_data_raw[LNneg_ERpos_HER2neg,]
expression_data=expression_data[,rownames(clinical_data)]

#Add in survival data - use DSS data (rather than OS data)
clinical_data[,"time_relapse"]=Complete_METABRIC_Clinical_Survival_Data__DSS[rownames(clinical_data),1]
clinical_data[,"time_relapse_yr"]=clinical_data[,"time_relapse"]/365
clinical_data[,"event_relapse"]=Complete_METABRIC_Clinical_Survival_Data__DSS[rownames(clinical_data),2]

#Create 10yr relapse variable similar to that used in training
time_relapse_yr=clinical_data[,"time_relapse_yr"]
event_relapse_10yr=clinical_data[,"event_relapse"]
event_relapse_10yr[which(time_relapse_yr<10 & event_relapse_10yr==0)]=NA
event_relapse_10yr[which(time_relapse_yr>10 & event_relapse_10yr==1)]=0
clinical_data[,"event_relapse_10yr"]=event_relapse_10yr

event_relapse_10yr_cens=clinical_data[,"event_relapse"]
event_relapse_10yr_cens[which(time_relapse_yr>10 & event_relapse_10yr_cens==1)]=0
clinical_data[,"event_relapse_10yr_cens"]=event_relapse_10yr_cens


#Extract only probes/genes needed for 17/8 gene models
#Need to choose representative probesets from Rosetta Chip
model_vars=rownames(rf_model$importance)
illumina_model_vars_all=names(probe_mapping[which(probe_mapping %in% model_vars)]) #All probes, including redundant for top8opt probes/genes

#probe_mapping[which(probe_mapping[,"GeneName"] %in% model_vars),c("RosettaProbeID","GeneName","EntrezID","AccNum")]
length(which(model_vars %in% probe_mapping))

#Use max variance to choose between probesets if more than one per gene
getmaxvarprobe=function(x){
 illumina_probes=as.vector(names(probe_mapping[which(probe_mapping==x)]))
 vars=apply(expression_data[illumina_probes,,drop=FALSE],1,var)
 maxvarprobe=names(which(vars==max(vars)))
 return(maxvarprobe)
}
illumina_model_vars=sapply(model_vars,getmaxvarprobe)

#Normalize/scale expression values to same distribution as used in training
#Determine target distribution based on cell line data
target_dist=normalize.quantiles.determine.target(as.matrix(training_data))

#Repeat above but using only probes/genes of interest
#top8opt_probes=c("10682_at","9133_at","1164_at","208079_s_at","23224_at","55435_at","23405_at","29127_at")
#top8opt_training_data=training_data[training_features[,"probes.ALL"]%in%top8opt_probes,]
#training_features_top8opt=training_features[training_features[,"probes.ALL"]%in%top8opt_probes,]
#rownames(top8opt_training_data)=training_features_top8opt[,"symbol.ALL"]
#target_dist_top8opt=normalize.quantiles.determine.target(as.matrix(top8opt_training_data))

#Set patient tumor data to target
expression_data_norm=normalize.quantiles.use.target(as.matrix(expression_data),target_dist)
rownames(expression_data_norm)=rownames(expression_data)
colnames(expression_data_norm)=colnames(expression_data)

#Repeat above using only probes/genes of interest
#expression_data_top8opt=expression_data[illumina_model_vars,]
#expression_data_top8opt_norm=normalize.quantiles.use.target(as.matrix(expression_data_top8opt),target_dist_top8opt)
#rownames(expression_data_top8opt_norm)=rownames(expression_data_top8opt)
#colnames(expression_data_top8opt_norm)=colnames(expression_data_top8opt)

#Another strategy would be to rescale each variable to the same range as those variables in training
expression_data_toscale=expression_data[illumina_model_vars_all,]
top17opt_probes=c("204767_s_at","10682_at","201291_s_at","9133_at","1164_at","208079_s_at","23224_at","55435_at","23220_at","201461_s_at","202709_at","57122_at","23405_at","201483_s_at","29127_at","204416_x_at","10628_at")
top17opt_training_data=training_data[training_features[,"probes.ALL"]%in%top17opt_probes,]
training_features_top17opt=training_features[training_features[,"probes.ALL"]%in%top17opt_probes,]
rownames(top17opt_training_data)=training_features_top17opt[,"symbol.ALL"]
top17opt_training_data_ranges=cbind(apply(top17opt_training_data,1,min),apply(top17opt_training_data,1,max))

expression_data_rescaled=expression_data_toscale
for (i in 1:length(illumina_model_vars_all)){
	illumina_id=illumina_model_vars_all[i]
	symbol=probe_mapping[illumina_id]
	expression_data_rescaled[illumina_id,]=rescale(expression_data_toscale[illumina_id,], newrange=range(top17opt_training_data_ranges[symbol,]))
}


#Build forest to assess redundant probes to identify best performing probe for each predictor gene
#predictor_data_assess=t(expression_data_norm[illumina_model_vars_all,])
predictor_data_assess=t(expression_data_rescaled[illumina_model_vars_all,]) ###USE ONLY WITH RESCALED OPTION BELOW###
colnames(predictor_data_assess)=illumina_model_vars_all
target=clinical_data[,"event_relapse"]
target[target==0]="NoRelapse"
target[target==1]="Relapse"
target=as.factor(target)
rf_model_assess=randomForest(x=predictor_data_assess, y=target, importance = TRUE, ntree = 10001, proximity=TRUE)

#Get importance measures
rf_importances=importance(rf_model_assess, scale=FALSE)
getmaximpprobe =function(x){
	illumina_probes=as.vector(names(probe_mapping[which(probe_mapping==x)]))
 	imps=rf_importances[illumina_probes,"MeanDecreaseGini"]
	maximpprobe=names(which(imps==max(imps)))
	if(length(illumina_probes)==1){maximpprobe=illumina_probes}
	return(maximpprobe)
}
illumina_model_vars_best=sapply(model_vars,getmaximpprobe)

#Set final predictor data for testing
#predictor_data=t(expression_data_norm[illumina_model_vars,])
#colnames(predictor_data)=names(illumina_model_vars)

#Try instead with data renormalized using only top8opt data
#predictor_data=t(expression_data_top8opt_norm[illumina_model_vars,])
#colnames(predictor_data)=names(illumina_model_vars)

#Use best probes where multiple to choose from
#predictor_data=t(expression_data_norm[illumina_model_vars_best,])
#colnames(predictor_data)=names(illumina_model_vars_best)

#Use best probes where multiple to choose from and rescaled data instead of renormalized
predictor_data=t(expression_data_rescaled[illumina_model_vars_best,])
colnames(predictor_data)=names(illumina_model_vars_best)


#Run test data through forest
RF_predictions_responses=predict(rf_model, predictor_data, type="response")
RF_predictions_probs=predict(rf_model, predictor_data, type="prob")
RF_predictions_votes=predict(rf_model, predictor_data, type="vote")

#Determine RF risk score according to previously determined thresholds
#low < 0.333
#0.333 >= int < 0.606 
#high >= 0.606

RF_risk_group=RF_predictions_votes[,"Relapse"]
RF_risk_group[RF_predictions_votes[,"Relapse"]<0.333]="Low"
RF_risk_group[RF_predictions_votes[,"Relapse"]>=0.333 & RF_predictions_votes[,"Relapse"]<0.606]="Int"
RF_risk_group[RF_predictions_votes[,"Relapse"]>=0.606]="High"

#Also break down into low, int, high simply according to same proportions as in training?
#Specifically, in training set there were 46.7% in Low, 38.6% in Int, and 14.9% in High risk groups
custom_lowN=round(length(RF_predictions_votes[,"Relapse"])*0.467)
custom_IntN=round(length(RF_predictions_votes[,"Relapse"])*0.386)
custom_HighN=length(RF_predictions_votes[,"Relapse"])-custom_lowN-custom_IntN

RF_risk_group_custom=RF_predictions_votes[,"Relapse"]
samples_sorted=names(sort(RF_predictions_votes[,"Relapse"]))
RF_risk_group_custom[samples_sorted[1:custom_lowN]]="Low"
RF_risk_group_custom[samples_sorted[(custom_lowN+1):(custom_lowN+custom_IntN)]]="Int"
RF_risk_group_custom[samples_sorted[(custom_lowN+custom_IntN+1):length(samples_sorted)]]="High"

#Yet another option is to use Mclust to choose new cutoffs as before
#Attempt to define cutoffs by mixed-model clustering of votes/probabilities
#Force mclust to break into 3 clusters 
x=as.numeric(RF_predictions_votes[,"Relapse"])
mclust_margins=Mclust(x, G=3)
summary(mclust_margins, x) #Gives you list of values returned
classification_margins=mclust_margins$classification
num_clusters=mclust_margins$G
table(classification_margins) #Shows not much better than using pre-defined cutoffs. Everything in intermediate group
#Choose cutoffs - use max and min of "middle" cluster (#2) to break into three components
mm_cutoff1=min(x[classification_margins==2])
mm_cutoff2=max(x[classification_margins==2])
RF_risk_group_MM=RF_predictions_votes[,"Relapse"]
RF_risk_group_MM[RF_predictions_votes[,"Relapse"]<mm_cutoff1]="Low"
RF_risk_group_MM[RF_predictions_votes[,"Relapse"]>=mm_cutoff1 & RF_predictions_votes[,"Relapse"]<mm_cutoff2]="Int"
RF_risk_group_MM[RF_predictions_votes[,"Relapse"]>=mm_cutoff2]="High"

#Join predictions with clinical data
clindata_plusRF=cbind(clinical_data,RF_predictions_responses,RF_predictions_votes,RF_risk_group_custom)
#clindata_plusRF=cbind(clinical_data,RF_predictions_responses,RF_predictions_votes,RF_risk_group)

#write results to file
write.table(clindata_plusRF,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Determine performance statistics
#Use only patients with 10 yr FU????
clindata_plusRF_10yrFU=clindata_plusRF[!is.na(clindata_plusRF[,"event_relapse_10yr"]),]

#Create ROC curve plot and calculate AUC
#Can use Relapse vote fractions fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling Response vs NoResponse

pred=prediction(clindata_plusRF[,"Relapse"],clindata_plusRF[,"event_relapse"])
#pred=prediction(clindata_plusRF_10yrFU[,"Relapse"],clindata_plusRF_10yrFU[,"event_relapse_10yr"])

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





#Create survival plot and statistics
#Calculate logrank survival statistic between groups
#Create new dataframe with just necessary data
surv_data=clindata_plusRF[,c("time_relapse_yr","event_relapse_10yr_cens","RF_risk_group_custom")]
#surv_data=clindata_plusRF[,c("time_relapse","event_relapse","RF_risk_group_custom")]
#surv_data=clindata_plusRF[,c("time_relapse","event_relapse","RF_risk_group")]

#create a survival object using data
surv_data.surv = with(surv_data, Surv(time_relapse_yr, event_relapse_10yr_cens==1))
#surv_data.surv = with(surv_data, Surv(time_relapse, event_relapse==1))
#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ RF_risk_group_custom, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

#Linear test p-value 
#Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
#See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
#recode  risk groups as 1,2,3
#surv_data_lin=clindata_plusRF[,c("time_relapse","event_relapse","RF_risk_group_custom")]
surv_data_lin=clindata_plusRF[,c("time_relapse_yr","event_relapse_10yr_cens","RF_risk_group_custom")]
surv_data_lin[,"RF_risk_group_custom"]=as.vector(surv_data_lin[,"RF_risk_group_custom"])
surv_data_lin[which(surv_data_lin[,"RF_risk_group_custom"]=="Low"),"RF_risk_group_custom"]=1
surv_data_lin[which(surv_data_lin[,"RF_risk_group_custom"]=="Int"),"RF_risk_group_custom"]=2
surv_data_lin[which(surv_data_lin[,"RF_risk_group_custom"]=="High"),"RF_risk_group_custom"]=3
surv_data_lin[,"RF_risk_group_custom"]=as.numeric(surv_data_lin[,"RF_risk_group_custom"])
survpvalue_linear=summary(coxph(Surv(time_relapse_yr, event_relapse_10yr_cens)~RF_risk_group_custom, data=surv_data_lin))$sctest[3]
#survpvalue_linear=summary(coxph(Surv(time_relapse, event_relapse)~RF_risk_group_custom, data=surv_data_lin))$sctest[3]
#survpvalue_linear = round(as.numeric(survpvalue_linear), digits=3)
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3, scientific=TRUE)

##Plot KM curve
krfit.by_RFgroup = survfit(surv_data.surv ~ RF_risk_group_custom, data = surv_data)
pdf(file=KMplotfile)
colors = rainbow(3)
title="Survival by RFRS - Not downsampled"
plot(krfit.by_RFgroup, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
#Set order of categories, categories are by default assigned colors alphabetically by survfit
groups=levels(surv_data[,"RF_risk_group_custom"]) #returns unique factor levels sorted alphabetically
names(colors)=groups
groups_custom=c("Low","Int","High")
colors_custom=colors[groups_custom]
group_sizes_custom=table(surv_data[,"RF_risk_group_custom"])[groups_custom]
#groups=unique(surv_data[,"RF_risk_group_custom"])
#group_sizes=table(surv_data[,"RF_risk_group_custom"])[groups]
legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue_linear, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors_custom,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()






#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
table(clindata_plusRF[,"event_relapse"])
#Keep all NoRelapses
NoRelapseCases=which(clindata_plusRF[,"event_relapse"]==0)
RelapseCases=which(clindata_plusRF[,"event_relapse"]==1)

#Downsample Relapse cases so that they represent only 15% of total cases (together with all non-10yr_relapse):  x / (237+x) = 0.15 [solving for x, ~42]
#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
#I=10 #for testing
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=8)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=42, replace = FALSE, prob = NULL)
 case_predictions_all_combined_down=clindata_plusRF[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group_custom"]=="Low","event_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group_custom"]=="Int","event_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_risk_group_custom"]=="High","event_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100

 #Calculate logrank survival statistic between groups
 #Create new dataframe with just necessary data
 surv_data=case_predictions_all_combined_down[,c("time_relapse","event_relapse","RF_risk_group_custom")]
 #create a survival object using data
 surv_data.surv = with(surv_data, Surv(time_relapse, event_relapse==1))
 #Calculate p-value
 survdifftest=survdiff(surv_data.surv ~ RF_risk_group_custom, data = surv_data)
 survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)

 #Linear test p-value 
 #Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
 #See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
 #recode  risk groups as 1,2,3
 surv_data[,"RF_risk_group_custom"]=as.vector(surv_data[,"RF_risk_group_custom"])
 surv_data[which(surv_data[,"RF_risk_group_custom"]=="Low"),"RF_risk_group_custom"]=1
 surv_data[which(surv_data[,"RF_risk_group_custom"]=="Int"),"RF_risk_group_custom"]=2
 surv_data[which(surv_data[,"RF_risk_group_custom"]=="High"),"RF_risk_group_custom"]=3
 surv_data[,"RF_risk_group_custom"]=as.numeric(surv_data[,"RF_risk_group_custom"])
 survpvalue_linear=summary(coxph(Surv(time_relapse, event_relapse)~RF_risk_group_custom, data=surv_data))$sctest[3]

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
write.table(downsampledata,file=case_pred_1000downsamp_KMpvalue_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
