library(ROCR)
#Set working directory and filenames for Input/output
setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/Oncotype/train_survival/")
#setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/Oncotype/test_survival/")
#setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/Oncotype/test_train_survival/")

#Read in data (expecting a tab-delimited file with header line and rownames)
datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/train_survival/standardCDF/ALL_gcrma.txt"
#datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_survival/standardCDF/ALL_gcrma.txt"
#datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_train_survival/standardCDF/ALL_gcrma.txt"

clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.train.2.txt"
#clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.test.2.txt"
#clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.txt"

#output files
ROC_pdffile="Oncotype_ROC.pdf"

combined_case_pred_downsamp_outfile="Cepheid_OncotypePredictions_combined_downsamp.txt"
combined_case_pred_downsamp_outfile2="Cepheid_OncotypePredictions_combined_downsamp.2.txt"

combined_case_pred_1000downsamp_10yrRFS_outfile="Cepheid_OncotypePredictions_combined_1000downsamp.txt"
combined_case_pred_1000downsamp_10yrRFS_outfile2="Cepheid_OncotypePredictions_combined_1000downsamp.2.txt"
oncotype_scores_outfile="Cepheid_OncotypePredictions_combined.txt"

data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

header=colnames(rawdata)

#Extract data need for Oncotype calculations
expr_data=rawdata[,4:length(header)]
probes=rawdata[,1]
rownames(expr_data)=probes
geneNames=as.vector(rawdata[,3])

###########################################################################################
#Calculate Oncotpye scores for data
##expected values are: Reference-normalized expression measurements range from 0 to 15, with a 1-unit increase reflecting approximately a doubling of RNA.
##this comes from the legend of figure 1 of the Paik et al., NEJM 2004
oncotypeDefinition <- function(){
 ##all genes with their groups
 genes <- list(her2=c('GRB7','ERBB2'),
               er=c('ESR1','PGR','BCL2','SCUBE2'),
               prolif=c('MKI67','STK15','BIRC5','CCNB1','MYBL2'),
               invas=c('MMP11','CTSL2'),
               cd68='CD68',gstm1='GSTM1',bag1='BAG1')
 controls <- c('ACTB','GAPDH','RPLPO','GUS','TFRC')
 ##weights of each gene within a group
 ##we will compute the mean of each afterward
 weights <- list(her2=c(0.9,0.1)*2,
                 er=c(0.8,1.2,1,1),
                 prolif=rep(1,5),
                 invas=rep(1,2),
                 cd68=1,gstm1=1,bag1=1)

 ##min weight of each group, 0 means no minimum
 mins <- c(er=0,her2=8,prolif=6.5, invas=0,cd68=0,gstm1=0,bag1=0)
 ##unscaled recurrence score, how the groups are added together
 rsu <- c(her=0.47, er=-0.34, prolif=1.04, invas=0.10,
          cd68=0.05, gstm1=-0.08, bag1=-0.07)
 ##RS=0 if RSU<0; RS=20×(RSU–6.7) if 0=RSU=100; and RS=100 if RSU>100.
 answer <- list(genes=genes,controls=controls,weights=weights,mins=mins,rsu=rsu)
}

oncotypeProbes <- function(geneNames,def=oncotypeDefinition(),
                          stat=apply(expr_data,1,mean)){
 ##get best probe(set) per sample
 genes.i <- sapply(def$genes,sapply,function(x){
   m <- which(geneNames%in%x)
   return(m[which.max(stat[m])])
 })
 controls.i <- sapply(def$controls,function(x){
   m <- which(geneNames%in%x)
   return(m[which.max(stat[m])])
 })
 def$genes.i <- genes.i
 def$controls.i <- controls.i
 def
}

##calculate oncotype for a given sample
computeOncotype <- function(expr,def=oncotypeProbes(geneNames),returnRaw=F){

 controls <- unlist(sapply(def$controls.i,sapply,function(x)mean(expr[x])))
 ##THIS IS COMPLETELY ARBITRARY (mean of controls over training set is 12.5)
 expr <- expr+16-mean(controls)

 ##list (groups) of values for each gene
 genes.raw <- sapply(def$genes.i,sapply,function(x)mean(expr[x]))
 ##apply weights
 genes.w<-sapply(1:length(genes.raw),function(x)genes.raw[[x]]*def$weights[[x]])
 ##we use the mean instead of a sum so that we can average out missing genes
 groups.raw <- sapply(genes.w,mean,na.rm=T)
 if (returnRaw)
   return(sum(groups.raw*def$rsu))
 groups.min <- pmax(groups.raw,def$mins)
 groups.w <- groups.min*def$rsu
 rsu <- 20*(sum(groups.w)-6.7)
 rsu
}

###########################################################################################

#Calculate Oncotype score for each sample
#test=expr_data[,1]
#names(test)=probes
#computeOncotype(test)
oncotype_scores=apply(expr_data,2,computeOncotype)

#Create ROC curve plot and calculate AUC
#Use Oncotype scores as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling recurrence vs non-recurrence
#First reduce to just patients with 10yr relapse
cases = which(!is.na(clindata[,"X10yr_relapse"]))
target=clindata[cases,"X10yr_relapse"]
predictions=as.vector(oncotype_scores[cases])
pred=prediction(predictions,target)
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


#Assign risk groups based on Oncotype score
#Low risk RS < 18, Intermediate risk RS >= 18 and < 31, High risk RS = 31
oncotype_group=oncotype_scores
oncotype_group[oncotype_scores<18]="low"
oncotype_group[oncotype_scores>=31]="high"
oncotype_group[oncotype_scores>=18 & oncotype_scores<31]="int"

#Rescale thresholds to produce similar % patients in each group to original Oncotype paper
#Determined by manual inspection
oncotype_group2=oncotype_scores
oncotype_group2[oncotype_scores<42]="low"
oncotype_group2[oncotype_scores>=56]="high"
oncotype_group2[oncotype_scores>=42 & oncotype_scores<56]="int"

#Add onto clinical data
clindata_oncotype=cbind(clindata,oncotype_scores,oncotype_group,oncotype_group2)
write.table(clindata_oncotype,file=oncotype_scores_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Keep all NoRelapses
NoRelapseCases=which(is.na(clindata_oncotype[,"X10yr_relapse"]) | clindata_oncotype[,"X10yr_relapse"]==0)
RelapseCases=which(clindata_oncotype[,"X10yr_relapse"]==1)

#Downsample Relapse cases so that they represent only 15% of total cases:  x / (645+x) = 0.15 [solving for x, ~114]
target_relapse_fraction=0.15
downsamp_relapse_target = round((length(NoRelapseCases)*target_relapse_fraction)/(1-target_relapse_fraction))

#Do multiple downsampling and get N and 10yr relapse rates for each risk group (oncotype_group), then average
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=6)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=downsamp_relapse_target, replace = FALSE, prob = NULL)
 case_predictions_all_combined_down=clindata_oncotype[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="low","X10yr_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="int","X10yr_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="high","X10yr_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
 downsampledata[i,1:6]=c(perc_10yr_relapse_low,length(low_10yr_relapses),perc_10yr_relapse_int,length(int_10yr_relapses),perc_10yr_relapse_high,length(high_10yr_relapses))
}
colnames(downsampledata)=c("low_perc","low_N","int_perc","int_N","high_perc","high_N")
#Print means to screen
mean(downsampledata[,"low_perc"]);mean(downsampledata[,"int_perc"]);mean(downsampledata[,"high_perc"])
write.table(downsampledata,file=combined_case_pred_1000downsamp_10yrRFS_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Do multiple downsampling and get N and 10yr relapse rates for each risk group (oncotype_group2), then average
I=1000
downsampledata2=matrix(data=NA, nrow=I, ncol=6)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=downsamp_relapse_target, replace = FALSE, prob = NULL)
 case_predictions_all_combined_down=clindata_oncotype[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="low","X10yr_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="int","X10yr_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="high","X10yr_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
 downsampledata2[i,1:6]=c(perc_10yr_relapse_low,length(low_10yr_relapses),perc_10yr_relapse_int,length(int_10yr_relapses),perc_10yr_relapse_high,length(high_10yr_relapses))
}
colnames(downsampledata2)=c("low_perc","low_N","int_perc","int_N","high_perc","high_N")
#Print means to screen
mean(downsampledata2[,"low_perc"]);mean(downsampledata2[,"int_perc"]);mean(downsampledata2[,"high_perc"])
write.table(downsampledata2,file=combined_case_pred_1000downsamp_10yrRFS_outfile2, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Create representative result for plotting survival figure (oncotype_group): run the block below until %'s similar to means above
random_RelapseCases=sample(x=RelapseCases, size=downsamp_relapse_target, replace = FALSE, prob = NULL)
case_predictions_all_combined_down=clindata_oncotype[c(NoRelapseCases,random_RelapseCases),]
low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="low","X10yr_relapse"]
int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="int","X10yr_relapse"]
high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group"]=="high","X10yr_relapse"]
perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
perc_10yr_relapse_low; perc_10yr_relapse_int; perc_10yr_relapse_high #Check values against means for above

write.table(case_predictions_all_combined_down,file=combined_case_pred_downsamp_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#Create representative result for plotting survival figure (oncotype_group2): run the block below until %'s similar to means above
random_RelapseCases=sample(x=RelapseCases, size=downsamp_relapse_target, replace = FALSE, prob = NULL)
case_predictions_all_combined_down=clindata_oncotype[c(NoRelapseCases,random_RelapseCases),]
low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="low","X10yr_relapse"]
int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="int","X10yr_relapse"]
high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"oncotype_group2"]=="high","X10yr_relapse"]
perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
perc_10yr_relapse_low; perc_10yr_relapse_int; perc_10yr_relapse_high #Check values against means for above

write.table(case_predictions_all_combined_down,file=combined_case_pred_downsamp_outfile2, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
