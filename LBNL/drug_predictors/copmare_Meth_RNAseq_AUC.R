#Compare AUC values for RNAseq (gene-level) vs methylation

AUC_datafile="/Users/ogriffit/Dropbox/drug_predictors/Rscripts/Rtoolbox/ExtraFiles/PredictorsMapping_PlatformAgnosticToolbox_RNAseq_AUC.txt"

AUCdata=read.table(AUC_datafile, header=TRUE, sep="\t", row.names=1)

cor.test(x=AUCdata[,1],y=AUCdata[,2],method="pearson")
cor.test(x=AUCdata[,1],y=AUCdata[,2],method="spearman")