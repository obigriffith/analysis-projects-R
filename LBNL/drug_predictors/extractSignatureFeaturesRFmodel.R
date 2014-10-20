#Load model to extract features from
RFmodeldir="/Users/ogriffit/Dropbox/drug_predictors/Rscripts/Rtoolbox/ModelsAgnostic/"
#RFmodelname="RF_model_RS_Bortezomib"
#RFmodelname="RF_model_RS_Cisplatin"
#RFmodelname="RF_model_EA_Meth_Ixabepilone"
#RFmodelname="RF_model_EA_GSK2119563A"
#RFmodelname="RF_model_U133A_SNP_GSK2 (PLKi)"
#RFmodelname="RF_model_EA_GSK2126458A"
#RFmodelname="RF_model_Meth_SNP_GSK2141795c"
#RFmodelname="RF_model_SNP_NU6102"
#RFmodelname="RF_model_EA_PF-3084014"
#RFmodelname="RF_model_U133A_Meth_SNP_Paclitaxel"
RFmodelname="RF_model_EA_AKT1-2 inhibitor"

RFmodelpath=paste(RFmodeldir,RFmodelname,".Rdata", sep="")

load(RFmodelpath)
gene_names=rownames(rf_model$importance)
#gene_names=gsub("__",",",gene_names)

sig_genes=cbind(gene_names,rf_model$importance[,"MeanDecreaseGini"],rf_model$importance[,"MeanDecreaseAccuracy"])
colnames(sig_genes)=c("Feature","VarImp","VarImpSD")


setwd("/Users/ogriffit/Dropbox/drug_predictors/TCGA/RtoolboxResults/Results/RFsig_features")
siggenesfilepath=paste(RFmodelname,".txt",sep="")
write.table(sig_genes, file=siggenesfilepath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)