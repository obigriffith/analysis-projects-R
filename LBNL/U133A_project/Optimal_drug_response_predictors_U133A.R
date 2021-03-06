#!/usr/bin/env Rscript

args=(commandArgs(TRUE))
CELfile=args[1]

library(Biobase)
library(affy)
library(gcrma)
# Custom cdf file (Dai et al, Michigan)
library(hgu133ahsentrezgcdf)
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)
library(randomForest)

path='/Users/adaemen/Dropbox/drug_predictors/U133A_project/'
setwd(path)
pathCELfiles='/Users/adaemen/Dropbox/drug_predictors/U133A_project/CEL_tumor/'

### Load tumor CEL file
txt_output=paste(CELfile,"_results.txt",sep="")
affy.data1=ReadAffy(filenames=paste(pathCELfiles,CELfile,sep=""))
setwd(paste(path,'CoreCELfiles/',sep=""))
affy.data2=ReadAffy()
# Combine tumor CEL file and cell line CEL files
affy.data=merge(affy.data1,affy.data2)

setwd(paste(path,'Results/',sep=""))

# Preprocessing: RMA, standard CDF
rma.standard.data=rma(affy.data)
exprSet.rma.standard=exprs(rma.standard.data)

# Preprocessing: GCRMA, standard CDF
gcrma.standard.data=gcrma(affy.data)
exprSet.gcrma.standard=exprs(gcrma.standard.data)

# Preprocessing: RMA, custom CDF
affy.data@cdfName="Hgu133A_Hs_ENTREZG"
rma.custom.data=rma(affy.data)
exprSet.rma.custom=exprs(rma.custom.data)
probes <- row.names(exprSet.rma.custom)
symbol <- unlist(mget(probes, hgu133ahsentrezgSYMBOL))
row.names(exprSet.rma.custom)=symbol

# Preprocessing: GCRMA, custom CDF
affy.data@cdfName="Hgu133A_Hs_ENTREZG"
gcrma.custom.data=gcrma(affy.data)
exprSet.gcrma.custom=exprs(gcrma.custom.data)
probes <- row.names(exprSet.gcrma.custom)
symbol <- unlist(mget(probes, hgu133ahsentrezgSYMBOL))
row.names(exprSet.gcrma.custom)=symbol


### Selection of stringent set of genes for each preprocessing
GeneProbeMapping.rma=read.table(paste(path,'ExtraFiles/Neve_AffyRMA_gene_maxvarprobe_mapping.txt',sep=""))
exprSet.rma.standard2=exprSet.rma.standard[as.vector(GeneProbeMapping.rma[-1,2]),]
rownames(exprSet.rma.standard2)=as.vector(GeneProbeMapping.rma[-1,1])
raw_data_import1=read.csv(paste(path,'ExtraFiles/Neve_AffyRMA_genelevel_maxvar_stringent.csv',sep=""), row.names=1)
genes.rma.standard=rownames(raw_data_import1)
expr.rma.standard=exprSet.rma.standard2[genes.rma.standard,1]

GeneProbeMapping.gcrma=read.table(paste(path,'ExtraFiles/Neve_AffyGCRMA_gene_maxvarprobe_mapping.txt',sep=""))
exprSet.gcrma.standard2=exprSet.gcrma.standard[as.vector(GeneProbeMapping.gcrma[-1,2]),]
rownames(exprSet.gcrma.standard2)=as.vector(GeneProbeMapping.gcrma[-1,1])
raw_data_import2=read.csv(paste(path,'ExtraFiles/Neve_AffyGCRMA_genelevel_maxvar_stringent.csv',sep=""), row.names=1)
genes.gcrma.standard=rownames(raw_data_import2)
expr.gcrma.standard=exprSet.gcrma.standard2[genes.gcrma.standard,1]

raw_data_import3=read.csv(paste(path,'ExtraFiles/Neve_AffyRMAcustom_genelevel_stringent.csv',sep=""), row.names=1)
genes.rma.custom=rownames(raw_data_import3)
expr.rma.custom=exprSet.rma.custom[genes.rma.custom,1]

raw_data_import4=read.csv(paste(path,'ExtraFiles/Neve_AffyGCRMAcustom_genelevel_stringent.csv',sep=""), row.names=1)
genes.gcrma.custom=rownames(raw_data_import4)
expr.gcrma.custom=exprSet.gcrma.custom[genes.gcrma.custom,1]


### Testing of best predictor per drug compound on tumor sample
drugs_interest=c("5-FU","AKT1-2 inhibitor","Cisplatin","Docetaxel","Vorinostat","Vinorelbine","GSK2126458A","GSK_Tykerb")
predictor_perf = data.frame(cbind(drugs_interest, Model=NA, Preprocessing=NA, Score=NA, Probability=NA, Sensitivity=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drugs_interest

# Cisplatin: LS-SVM - GCRMA standard CDF
drug="Cisplatin"
ModelCisplatin=read.table(paste(path,'Models/Model_GCRMA_standardCDF_Cisplatin.txt',sep=""),as.is=1,header=TRUE)
TrainCisplatin=read.table(paste(path,'Models/Trainingdata_GCRMA_standardCDF_Cisplatin.txt',sep=""),header=TRUE)
Asigm=ModelCisplatin[[1]][1]
Bsigm=ModelCisplatin[[1]][2]
bcoeff=ModelCisplatin[[1]][3]
alpha=ModelCisplatin[[1]][-(1:3)]
target=as.matrix(TrainCisplatin[1,])
traindata=as.matrix(TrainCisplatin[-1,])

Ktest=t(traindata) %*% as.matrix(expr.gcrma.standard)
latent=((alpha*target) %*% Ktest) + bcoeff
predictor_perf[drug,"Score"]=format(latent,digits=5)
prob=1/(1+exp(Asigm*latent+Bsigm))
predictor_perf[drug,"Probability"]=format(prob,digits=5)
predictor_perf[drug,"Sensitivity"]=(prob>0.5)
predictor_perf[drug,"Model"]="LS-SVM"
predictor_perf[drug,"Preprocessing"]="GCRMA standard CDF"

# Docetaxel: LS-SVM - GCRMA custom CDF
drug="Docetaxel"
ModelDocetaxel=read.table(paste(path,'Models/Model_GCRMA_customCDF_Docetaxel.txt',sep=""),as.is=1,header=TRUE)
TrainDocetaxel=read.table(paste(path,'Models/Trainingdata_GCRMA_customCDF_Docetaxel.txt',sep=""),header=TRUE)
Asigm=ModelDocetaxel[[1]][1]
Bsigm=ModelDocetaxel[[1]][2]
bcoeff=ModelDocetaxel[[1]][3]
alpha=ModelDocetaxel[[1]][-(1:3)]
target=as.matrix(TrainDocetaxel[1,])
traindata=as.matrix(TrainDocetaxel[-1,])

Ktest=t(traindata) %*% as.matrix(expr.gcrma.custom)
latent=((alpha*target) %*% Ktest) + bcoeff
predictor_perf[drug,"Score"]=format(latent,digits=5)
prob=1/(1+exp(Asigm*latent+Bsigm))
predictor_perf[drug,"Probability"]=format(prob,digits=5)
predictor_perf[drug,"Sensitivity"]=(prob>0.5)
predictor_perf[drug,"Model"]="LS-SVM"
predictor_perf[drug,"Preprocessing"]="GCRMA custom CDF"

# GSK2126458A: LS-SVM - RMA custom CDF
drug="GSK2126458A"
ModelGSK212=read.table(paste(path,'Models/Model_RMA_customCDF_GSK2126458A.txt',sep=""),as.is=1,header=TRUE)
TrainGSK212=read.table(paste(path,'Models/Trainingdata_RMA_customCDF_GSK2126458A.txt',sep=""),header=TRUE)
Asigm=ModelGSK212[[1]][1]
Bsigm=ModelGSK212[[1]][2]
bcoeff=ModelGSK212[[1]][3]
alpha=ModelGSK212[[1]][-(1:3)]
target=as.matrix(TrainGSK212[1,])
traindata=as.matrix(TrainGSK212[-1,])

Ktest=t(traindata) %*% as.matrix(expr.rma.custom)
latent=((alpha*target) %*% Ktest) + bcoeff
predictor_perf[drug,"Score"]=format(latent,digits=5)
prob=1/(1+exp(Asigm*latent+Bsigm))
predictor_perf[drug,"Probability"]=format(prob,digits=5)
predictor_perf[drug,"Sensitivity"]=(prob>0.5)
predictor_perf[drug,"Model"]="LS-SVM"
predictor_perf[drug,"Preprocessing"]="RMA custom CDF"

# GSK_Tykerb: LS-SVM - RMA custom CDF
drug="GSK_Tykerb"
ModelGSKTykerb=read.table(paste(path,'Models/Model_RMA_customCDF_GSK_Tykerb.txt',sep=""),as.is=1,header=TRUE)
TrainGSKTykerb=read.table(paste(path,'Models/Trainingdata_RMA_customCDF_GSK_Tykerb.txt',sep=""),header=TRUE)
Asigm=ModelGSKTykerb[[1]][1]
Bsigm=ModelGSKTykerb[[1]][2]
bcoeff=ModelGSKTykerb[[1]][3]
alpha=ModelGSKTykerb[[1]][-(1:3)]
target=as.matrix(TrainGSKTykerb[1,])
traindata=as.matrix(TrainGSKTykerb[-1,])

Ktest=t(traindata) %*% as.matrix(expr.rma.custom)
latent=((alpha*target) %*% Ktest) + bcoeff
predictor_perf[drug,"Score"]=format(latent,digits=5)
prob=1/(1+exp(Asigm*latent+Bsigm))
predictor_perf[drug,"Probability"]=format(prob,digits=5)
predictor_perf[drug,"Sensitivity"]=(prob>0.5)
predictor_perf[drug,"Model"]="LS-SVM"
predictor_perf[drug,"Preprocessing"]="RMA custom CDF"

# AKT1-2 inhibitor: RF - RMA standard CDF
drug="AKT1-2 inhibitor"
RFdata=as.data.frame(expr.rma.standard)
colnames(RFdata)=CELfile
RFdata=t(RFdata)
load(file=paste(path,'Models/RFModel_RMA_standardCDF_AKT12inhibitor',sep=""))
RF_predictions_response=predict(rf_model, RFdata, type="response")
RF_predictions_prob=predict(rf_model, RFdata, type="prob")
RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
predictor_perf[drug,"Score"]=RF_predictions_vote[1,"sensitive"]
predictor_perf[drug,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
predictor_perf[drug,"Sensitivity"]=RF_predictions_response=="sensitive"
predictor_perf[drug,"Model"]="RF"
predictor_perf[drug,"Preprocessing"]="RMA standard CDF"

# Vinorelbine: RF - GCRMA custom CDF
drug="Vinorelbine"
RFdata=as.data.frame(expr.gcrma.custom)
colnames(RFdata)=CELfile
RFdata=t(RFdata)
load(file=paste(path,'Models/RFModel_GCRMA_customCDF_Vinorelbine',sep=""))
RF_predictions_response=predict(rf_model, RFdata, type="response")
RF_predictions_prob=predict(rf_model, RFdata, type="prob")
RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
predictor_perf[drug,"Score"]=RF_predictions_vote[1,"sensitive"]
predictor_perf[drug,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
predictor_perf[drug,"Sensitivity"]=RF_predictions_response=="sensitive"
predictor_perf[drug,"Model"]="RF"
predictor_perf[drug,"Preprocessing"]="GCRMA custom CDF"

# Vorinostat: RF - GCRMA custom CDF
drug="Vorinostat"
RFdata=as.data.frame(expr.gcrma.custom)
colnames(RFdata)=CELfile
RFdata=t(RFdata)
load(file=paste(path,'Models/RFModel_GCRMA_customCDF_Vorinostat',sep=""))
RF_predictions_response=predict(rf_model, RFdata, type="response")
RF_predictions_prob=predict(rf_model, RFdata, type="prob")
RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
predictor_perf[drug,"Score"]=RF_predictions_vote[1,"sensitive"]
predictor_perf[drug,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
predictor_perf[drug,"Sensitivity"]=RF_predictions_response=="sensitive"
predictor_perf[drug,"Model"]="RF"
predictor_perf[drug,"Preprocessing"]="GCRMA custom CDF"

# 5-FU: RF - GCRMA standard CDF
drug="5-FU"
RFdata=as.data.frame(expr.gcrma.standard)
colnames(RFdata)=CELfile
RFdata=t(RFdata)
load(file=paste(path,'Models/RFModel_GCRMA_standardCDF_5FU',sep=""))
RF_predictions_response=predict(rf_model, RFdata, type="response")
RF_predictions_prob=predict(rf_model, RFdata, type="prob")
RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
predictor_perf[drug,"Score"]=RF_predictions_vote[1,"sensitive"]
predictor_perf[drug,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
predictor_perf[drug,"Sensitivity"]=RF_predictions_response=="sensitive"
predictor_perf[drug,"Model"]="RF"
predictor_perf[drug,"Preprocessing"]="GCRMA standard CDF"


### Save results to txt file
write.table(predictor_perf, file=txt_output, sep="\t", row.names=FALSE)