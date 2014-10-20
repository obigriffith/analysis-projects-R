library(ROCR)
library(dplyr)
#Set working directory and filenames for Input/output
setwd("~/garibaldi/meissto/repos/brcarecurrence")
#setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/Oncotype/test_survival/")
#setwd("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/Oncotype/test_train_survival/")

#Read in data (expecting a tab-delimited file with header line and rownames)
datafile_train1="data/processed_final2_train_survival_combined_gcrma.1.txt"
datafile_train2="data/processed_final2_train_survival_combined_gcrma.2.txt"
datafile_test="data/processed_final2_test_survival_combined_gcrma.txt"
#datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_survival/standardCDF/ALL_gcrma.txt"
#datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_train_survival/standardCDF/ALL_gcrma.txt"

clindatafile_train="data/clinical_training.txt"
clindatafile_test="data/clinical_test.txt"
#clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.test.2.txt"
#clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.txt"

#output files
ROC_pdffile="Oncotype_ROC.pdf"

combined_case_pred_downsamp_outfile="Cepheid_OncotypePredictions_combined_downsamp.txt"
combined_case_pred_downsamp_outfile2="Cepheid_OncotypePredictions_combined_downsamp.2.txt"

combined_case_pred_1000downsamp_10yrRFS_outfile="Cepheid_OncotypePredictions_combined_1000downsamp.txt"
combined_case_pred_1000downsamp_10yrRFS_outfile2="Cepheid_OncotypePredictions_combined_1000downsamp.2.txt"
oncotype_scores_outfile="Cepheid_OncotypePredictions_combined.txt"

data_import_train1=read.table(datafile_train1, header = TRUE, na.strings = "NA", sep="\t")
data_import_train2=read.table(datafile_train2, header = FALSE, na.strings = "NA", sep="\t")
colnames(data_import_train2) <- colnames(data_import_train1)
data_import_train=rbind(data_import_train1, data_import_train2)

data_import_test=read.table(datafile_test, header = TRUE, na.strings = "NA", sep="\t")

data_import <- cbind(data_import_train, data_import_test[,-c(1:3)])

# remove duplicated probes on the u133b chip (2nd appearance)
x <- data_import$probes.ALL[which(duplicated(data_import$probes.ALL))]
rem <- which(data_import$probes.ALL %in% x)[-c(1:length(which(duplicated(data_import$probes.ALL))))]
data_import <- data_import[-rem, ]

clin_data_import_train=read.table(clindatafile_train, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import_test=read.table(clindatafile_test, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import <- rbind(clin_data_import_train, clin_data_import_test)

# THIS SORTING SEEMS NOT TO WORK!!!
# new code line 67-70
# ----------------
#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
#clin_data_order=order(clin_data_import[,"GSM"])
#clindata=clin_data_import[clin_data_order,]
#data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
#rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
# ----------------

#Extract data need for Oncotype calculations
expr_data=data_import[,4:dim(data_import)[2]]
probes=data_import[,1]
rownames(expr_data)=probes
geneNames=as.vector(data_import[,3])

# sort the clinical data
x <- colnames(expr_data)
x <- gsub('.CEL', '', x)
x <- gsub('.cel', '', x)
clindata <- clin_data_import[match(x, clin_data_import$GSM),]

###########################################################################################
#Calculate Oncotpye scores for data
##expected values are: Reference-normalized expression measurements range from 0 to 15, with a 1-unit increase reflecting approximately a doubling of RNA.
##this comes from the legend of figure 1 of the Paik et al., NEJM 2004

###############################################################################
# OncotypeDX function adapted from genefun package
###############################################################################

# function to compute oncotypeDX score
oncotypeDX <- function(x) {
  
  geneids <- c('4288', '6790', '332', '891', '4605', '2099', '5241', '596', '57758',
               '2886', '2064', '4320', '1515', '968', '2944', '573')
  
  ref <- mean(x[c('60', '2597', '6175', '2990', '7037')])
  
  norm <- x+16-ref # this seems to work better then the scaling from 0-15
  # as suggested in the original publication
  
  #norm <- scale(x[geneids], center=FALSE, scale=ref)
  
  # normScaled <- apply(norm, 2, function(x) {
  #   xx <- (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
  #                                       min(x, na.rm = TRUE))
  #   return(xx * 15)
  # })
  
  normScaled <- as.matrix(norm)
  
  grb7Group <- as.numeric(0.9 * normScaled['2886',] + 0.1 * normScaled['2064',])
  if (grb7Group < 8) grb7Group <- 8
  
  erGroup <- as.numeric((0.8 * normScaled['2099',] + 1.2 * normScaled['5241',] + normScaled['596',] + normScaled['57758',]) / 4)
  
  prolifGroup <- as.numeric((normScaled['332',] + normScaled['4288',] + normScaled['4605',] + normScaled['891',] + normScaled['6790',]) / 5)
  if (prolifGroup < 6.5) prolifGroup <- 6.5 
  
  invGroup <- as.numeric((normScaled['1515',] + normScaled['4320',]) / 2)
  
  cd68 <- as.numeric(normScaled['968', ])
  gstm1 <- as.numeric(normScaled['2944', ])
  bag1 <- as.numeric(normScaled['573', ])
  
  # unscaled recurance score
  RSu <- (0.47 * grb7Group) - (0.34 * erGroup) + (1.04 * prolifGroup) + (0.1 * invGroup) + (0.05 * cd68) - (0.08 * gstm1) - (0.07 * bag1)
  
  # scaled recucance score
  if (RSu < 0) {
    RS <- 0
  }
  else if (RSu >= 0 & RSu <= 100) {
    RS <- 20 * (RSu - 6.7)
  }
  else if (RSu > 100) {
    RS <- 100
  }
  
  return(list(RSu=RSu, RS=RS, grb7Group=grb7Group, erGroup=erGroup, prolifGroup=prolifGroup, invGroup=invGroup, cd68=cd68, gstm1=gstm1, bag1=bag1))
}

###########################################################################################

genes <- c('MKI67', 'AURKA', 'BIRC5', 'CCNB1', 'MYBL2', 'ESR1', 'PGR', 'BCL2', 'SCUBE2', 'GRB7',
           'ERBB2', 'MMP11', 'CTSL2', 'CD68', 'GSTM1', 'BAG1', 'ACTB', 'GAPDH', 'RPLP0', 'GUSB', 
           'TFRC')
geneids <- c('4288', '6790', '332', '891', '4605', '2099', '5241', '596', '57758',
             '2886', '2064', '4320', '1515', '968', '2944', '573','60', '2597', 
             '6175', '2990', '7037')

# find optimal probesets
# def. optimal probeset: highest mean exppression across all samples
stat=apply(expr_data,1,mean)
optProbes <- names(unlist(lapply(genes, function(x) which.max(stat[which(geneNames%in%x)]))))
optIn <- expr_data[optProbes, ]
rownames(optIn) <- geneids

#Calculate Oncotype score for each sample
oncodx <- apply(optIn,2,oncotypeDX)
oncotype_scores <- unlist(oncodx)[seq(2, length(unlist(oncodx)), by=9)]

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

# note: combining test+validation data, values change from the original code implemented on
# the test data only
oncotype_group2=oncotype_scores
oncotype_group2[oncotype_scores<57]="low"
oncotype_group2[oncotype_scores>=77]="high"
oncotype_group2[oncotype_scores>=57 & oncotype_scores<77]="int"

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


## KM Plot
library(survival)
clindata$odx <-oncotype_group2

fit <- survfit(Surv(t_rfs,e_rfs)~odx, data=clindata)
survdiff(Surv(t_rfs,e_rfs)~odx, data=clindata)
plot(fit)

## add group scores to clinical data
clindata$RSu <- unlist(oncodx)[seq(1, length(unlist(oncodx)), by=9)]
clindata$RS <- unlist(oncodx)[seq(2, length(unlist(oncodx)), by=9)]
clindata$grb7Group <- unlist(oncodx)[seq(3, length(unlist(oncodx)), by=9)]
clindata$erGroup <- unlist(oncodx)[seq(4, length(unlist(oncodx)), by=9)]
clindata$prolifGroup <- unlist(oncodx)[seq(5, length(unlist(oncodx)), by=9)]
clindata$invGroup <- unlist(oncodx)[seq(6, length(unlist(oncodx)), by=9)]
clindata$cd68 <- unlist(oncodx)[seq(7, length(unlist(oncodx)), by=9)]
clindata$gstm1 <- unlist(oncodx)[seq(8, length(unlist(oncodx)), by=9)]
clindata$bag1 <- unlist(oncodx)[seq(9, length(unlist(oncodx)), by=9)]

write.table(clindata,file='clindata_oncotype.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# ----------------------------------------------------------------------------
# prepare data for upload into branch
# only one probeset per gene (entrez_id)
# merge clin + expr data, keep only samples with 10y relapse info

allGenes <- unique(na.omit(data_import$ID.ALL))
optProbesAllGenes <- names(unlist(lapply(allGenes, function(x) which.max(stat[which(data_import$ID.ALL%in%x)]))))
expr_branch <- expr_data[optProbesAllGenes, ]
colnames(expr_branch) <- gsub('.CEL', '', colnames(expr_branch))
colnames(expr_branch) <- gsub('.cel', '', colnames(expr_branch))

branch_data_griffith <- cbind(clindata, t(expr_branch), X10yrsurvival=ifelse(clindata$X10yr_relapse==1, 'yes', 'no'))
branch_data_griffith <- branch_data_griffith[which(branch_data_griffith$X10yrsurvival=='yes' | branch_data_griffith$X10yrsurvival=='no'), ]

write.table(branch_data_griffith,file='branch_griffith.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#split up test/validation
branch_data_griffith_train <- branch_data_griffith[na.omit(match(clin_data_import_train$GSM, branch_data_griffith$GSM)),]
branch_data_griffith_test <- branch_data_griffith[na.omit(match(clin_data_import_test$GSM, branch_data_griffith$GSM)),]

write.table(branch_data_griffith_train,file='branch_griffith_train.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(branch_data_griffith_test,file='branch_griffith_test.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#mapping file entrez <-> symbol
mapping <- cbind(entrezid=allGenes, 
                 symbol=as.vector(data_import$symbol.ALL)[match(allGenes, data_import$ID.ALL)], 
                 probeset=optProbesAllGenes)
write.table(mapping,file='branch_griffith_mapping.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#-------------------------------------------------------------------------------------------
### Test on METABRIC data
#-------------------------------------------------------------------------------------------
test <- read.csv2('data/Metabric_clinical_expression_DSS_sample_filtered.tsv', sep='\t')
clinicalMeta <- read.csv2('data/Metabric_full_clinical_feats_sample_filtered.tsv', sep='\t')
survivalMeta <- read.csv2('data/Metabric_DS_survival_filtered.tsv', sep='\t')

exprMeta <- t(as.matrix(test[,14:26030]))
mode(exprMeta) <- 'numeric'
colnames(exprMeta) <- test$id
rm(test) # free up some memory..

#library(preprocessCore)
library(illuminaHumanv4.db)
#exprMeta.norm <- normalize.quantiles(exprMeta)
#dimnames(exprMeta.norm) <- dimnames(exprMeta)

ill_entrez <- unlist(mget(rownames(exprMeta), illuminaHumanv4ENTREZID))

bestIll <- function(gid) {
  probes <- unlist(mget(gid, illuminaHumanv4ALIAS2PROBE))
  mat <- exprMeta[which(rownames(exprMeta) %in% probes), ,drop=FALSE]
  mat_mean <- apply(mat, 1, mean)
  mat_best <- mat[which.max(mat_mean), , drop=F]
  return(mat_best)
}

# note: BAG1 seems to be HAP (older symbol) on illumina annotation
genes <- c('MKI67', 'AURKA', 'BIRC5', 'CCNB1', 'MYBL2', 'ESR1', 'PGR', 'BCL2', 'SCUBE2', 'GRB7',
           'ERBB2', 'MMP11', 'CTSV', 'CD68', 'GSTM1', 'HAP', 'ACTB', 'GAPDH', 'RPLP0', 'GUSB', 
           'TFRC')
xx <- lapply(genes, bestIll)
meta_in <- do.call(rbind, xx)
ill_odx_probests <- rownames(meta_in)
rownames(meta_in) <- geneids

metaDX <- apply(meta_in, 2, oncotypeDX)
metaDXRS <- unlist(metaDX)[seq(2, length(unlist(metaDX)), by=9)]

# '10 y relapse'
survivalMeta$y10 <- ifelse(survivalMeta$time/365.25>10, 0, 1)

predictions=as.vector(metaDXRS)
pred=prediction(predictions,survivalMeta$y10)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
AUC_out=paste("AUC=",AUC,sep="")
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))

metaDXrs_new=metaDXRS
metaDXrs_new[metaDXrs_new<131]=0
metaDXrs_new[metaDXrs_new>=131 & metaDXrs_new<145]=0.5
metaDXrs_new[metaDXrs_new>=145]=1

survivalMeta$odx <-metaDXrs_new

fit <- survfit(Surv(time,status)~odx, data=survivalMeta)
survdiff(Surv(time,status)~odx, data=survivalMeta)
plot(fit)

## add group scores to clinical data
survivalMeta$RSu <- unlist(metaDX)[seq(1, length(unlist(metaDX)), by=9)]
survivalMeta$RS <- unlist(metaDX)[seq(2, length(unlist(metaDX)), by=9)]
survivalMeta$grb7Group <- unlist(metaDX)[seq(3, length(unlist(metaDX)), by=9)]
survivalMeta$erGroup <- unlist(metaDX)[seq(4, length(unlist(metaDX)), by=9)]
survivalMeta$prolifGroup <- unlist(metaDX)[seq(5, length(unlist(metaDX)), by=9)]
survivalMeta$invGroup <- unlist(metaDX)[seq(6, length(unlist(metaDX)), by=9)]
survivalMeta$cd68 <- unlist(metaDX)[seq(7, length(unlist(metaDX)), by=9)]
survivalMeta$gstm1 <- unlist(metaDX)[seq(8, length(unlist(metaDX)), by=9)]
survivalMeta$bag1 <- unlist(metaDX)[seq(9, length(unlist(metaDX)), by=9)]

survivalMeta$sample <- rownames(survivalMeta)
write.table(survivalMeta,file='oncotype_metabric.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# ----------------------------------------------------------------------------
# prepare data for upload into branch
# only one probeset per gene
# merge clin + expr data, keep only samples with 10y relapse info

# merge clinical + surv data
mergeClinMeta <- merge(clinicalMeta, survivalMeta, by.x='X', by.y='sample')

allGenes <- as.vector(unique(na.omit(ill_entrez)))
stat <- apply(exprMeta,1,mean)
optProbesAllGenes <- names(unlist(lapply(allGenes, function(x) which.max(stat[which(ill_entrez%in%x)]))))
expr_branch_metabric <- exprMeta[optProbesAllGenes, ]

branch_data_metabric <- cbind(mergeClinMeta, t(expr_branch_metabric), X10yrsurvival=ifelse(mergeClinMeta$y10==1, 'yes', 'no'))

write.table(branch_data_metabric,file='branch_metabric.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#mapping file entrez <-> symbol
mapping <- cbind(entrezid=allGenes, 
                 symbol=as.vector(unlist(mget(optProbesAllGenes, illuminaHumanv4SYMBOL))),
                 probeset=optProbesAllGenes)
write.table(mapping,file='branch_metabric_mapping.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#-------------------------------------------------------------------------------------------
### Test on Oslo data
#-------------------------------------------------------------------------------------------
oslo <- read.csv2('data/Oslo_clinical_expression_OS_sample_filt.tsv', sep='\t')
survivalOslo <- read.csv2('data/Oslo_survival_sample_filt.tsv', sep='\t')
clinicalOslo <- oslo[,1:13]

exprOslo <- t(as.matrix(oslo[,14:26030]))
mode(exprOslo) <- 'numeric'
colnames(exprMeta) <- oslo$id
rm(oslo) # free up some memory..

oslo_in <- exprOslo[ill_odx_probests,]
rownames(oslo_in) <-geneids

osloDX <- apply(oslo_in, 2, oncotypeDX)
osloDXRS <- unlist(osloDX)[seq(2, length(unlist(osloDX)), by=9)]

# '10 y relapse'
survivalOslo$y10 <- ifelse(survivalOslo$time/365.25>10, 0, 1)

predictions=as.vector(osloDXRS)
pred=prediction(predictions,survivalOslo$y10)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
AUC_out=paste("AUC=",AUC,sep="")
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))

osloDXrs_new=osloDXRS
osloDXrs_new[osloDXrs_new<131]=0
osloDXrs_new[osloDXrs_new>=131 & osloDXrs_new<145]=0.5
osloDXrs_new[osloDXrs_new>=145]=1

survivalOslo$odx <-osloDXrs_new

fit <- survfit(Surv(time,status)~odx, data=survivalOslo)
survdiff(Surv(time,status)~odx, data=survivalOslo)
plot(fit)

## add group scores to clinical data
survivalOslo$RSu <- unlist(osloDX)[seq(1, length(unlist(osloDX)), by=9)]
survivalOslo$RS <- unlist(osloDX)[seq(2, length(unlist(osloDX)), by=9)]
survivalOslo$grb7Group <- unlist(osloDX)[seq(3, length(unlist(osloDX)), by=9)]
survivalOslo$erGroup <- unlist(osloDX)[seq(4, length(unlist(osloDX)), by=9)]
survivalOslo$prolifGroup <- unlist(osloDX)[seq(5, length(unlist(osloDX)), by=9)]
survivalOslo$invGroup <- unlist(osloDX)[seq(6, length(unlist(osloDX)), by=9)]
survivalOslo$cd68 <- unlist(osloDX)[seq(7, length(unlist(osloDX)), by=9)]
survivalOslo$gstm1 <- unlist(osloDX)[seq(8, length(unlist(osloDX)), by=9)]
survivalOslo$bag1 <- unlist(osloDX)[seq(9, length(unlist(osloDX)), by=9)]

survivalOslo$sample <- rownames(survivalOslo)
write.table(survivalOslo,file='oncotype_oslo.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# ----------------------------------------------------------------------------
# prepare data for upload into branch
# only one probeset per gene
# merge clin + expr data, keep only samples with 10y relapse info

# merge clinical + surv data
mergeClinOslo <- merge(clinicalOslo, survivalOslo, by.x='id', by.y='sample')

expr_branch_oslo <- exprOslo[optProbesAllGenes, ]

branch_data_oslo <- cbind(mergeClinOslo, t(expr_branch_oslo), X10yrsurvival=ifelse(mergeClinOslo$y10==1, 'yes', 'no'))

write.table(branch_data_oslo,file='branch_oslo.csv', sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
