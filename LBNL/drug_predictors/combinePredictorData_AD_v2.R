#Combine data files for predictor analysis

outdir="/Users/adaemen/Dropbox/drug_predictors/combined"
setwd(outdir)
outfile="BCCL_combined_data.4.txt"
outfile2="BCCL_combined_data.4.noNA.txt"

#Cell line info
cell_line_datafile="/Users/adaemen/Dropbox/drug_predictors/BCCL_data_list_v2.txt"
cell_data_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=1)
cell_lines=sort(cell_data_import[,1])

### RNAseq data ###
RNAdatafile="/Users/adaemen/Dropbox/drug_predictors/RNAseq/filtered_v3/breastRNAseq_all_stringent.txt"
RNA_annotation_file="/Users/adaemen/Dropbox/drug_predictors/RNAseq/hs_53_36o_ENSG2Symbol.txt"
RNA_data_import=read.table(RNAdatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
RNA_anno_import=read.table(RNA_annotation_file, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2), row.names=1)
RNA_feat_data=RNA_data_import[,1:4]

#Look up Symbol for ENSG ID and add to feature data
ENSG2Symbol=RNA_anno_import[RNA_feat_data[,4],]
RNA_feat_data=cbind(RNA_feat_data,ENSG2Symbol)

#Create single combined ID for RNAseq: FID__SeqName
comb_ID=paste(RNA_feat_data[,2],RNA_feat_data[,3],RNA_feat_data[,5],sep="__")
RNA_var_names=cbind(rep("RS",length(comb_ID)),comb_ID)
colnames(RNA_var_names)=c("DataType","ID")
RNA_data=RNA_data_import[,5:length(colnames(RNA_data_import))]
RNA_data=log2(RNA_data+1) #Log2 data

#Fix misnamed libraries
colnames(RNA_data)[which(colnames(RNA_data)=="X184A1")]="184A1"
colnames(RNA_data)[which(colnames(RNA_data)=="X184B5")]="184B5"
colnames(RNA_data)[which(colnames(RNA_data)=="X600MPE")]="600MPE"
colnames(RNA_data)[which(colnames(RNA_data)=="X21MT1")]="21MT1"
colnames(RNA_data)[which(colnames(RNA_data)=="X21MT2")]="21MT2"
colnames(RNA_data)[which(colnames(RNA_data)=="X21NT")]="21NT"
colnames(RNA_data)[which(colnames(RNA_data)=="X21PT")]="21PT"

#Create NA columns for cell lines not represented
RNA_NA_lines=cell_lines[which(!cell_lines %in% colnames(RNA_data))]
RNA_NA_data=matrix(NA,nrow=length(comb_ID), ncol=length(RNA_NA_lines))
colnames(RNA_NA_data)=RNA_NA_lines
RNA_data_wNA=cbind(RNA_data,RNA_NA_data)[,cell_lines]
RNA_data_complete=cbind(RNA_var_names,RNA_data_wNA)


### exonarray data ###
EAdatafile="/Users/adaemen/Dropbox/drug_predictors/ExonArray/filtered_v2/breastExon_all_extrastringent.csv"
EA_data_import=read.csv(EAdatafile) 
EA_feat_data=as.vector(EA_data_import[,1])
EA_var_names=cbind(rep("EA",length(EA_feat_data)),EA_feat_data)
colnames(EA_var_names)=c("DataType","ID")
EA_data=EA_data_import[,2:length(colnames(EA_data_import))]

#Fix misnamed libraries
colnames(EA_data)[which(colnames(EA_data)=="X184A1")]="184A1"
colnames(EA_data)[which(colnames(EA_data)=="X184B5")]="184B5"
colnames(EA_data)[which(colnames(EA_data)=="X600MPE")]="600MPE"

#Create NA columns for cell lines not represented
EA_NA_lines=cell_lines[which(!cell_lines %in% colnames(EA_data))]
EA_NA_data=matrix(NA,nrow=length(EA_feat_data), ncol=length(EA_NA_lines))
colnames(EA_NA_data)=EA_NA_lines
EA_data_wNA=cbind(EA_data,EA_NA_data)[,cell_lines]
EA_data_complete=cbind(EA_var_names,EA_data_wNA)


### exome-seq mutation data ###
#ESdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/mutations/Subset5_filt_mutations.txt"
#ES_data_import=read.table(ESdatafile, header = TRUE, na.strings = "NA", sep="\t")
#ES_feat_data=as.vector(ES_data_import[,1])
#ES_var_names=cbind(rep("ES",length(ES_feat_data)),ES_feat_data)
#colnames(ES_var_names)=c("DataType","ID")
#ES_data=ES_data_import[,2:length(colnames(ES_data_import))]

#Fix misnamed libraries
#colnames(ES_data)[which(colnames(ES_data)=="SUM1315")]="SUM1315MO2"

#Create NA columns for cell lines not represented
#ES_NA_lines=cell_lines[which(!cell_lines %in% colnames(ES_data))]
#ES_NA_data=matrix(NA,nrow=length(ES_feat_data), ncol=length(ES_NA_lines))
#colnames(ES_NA_data)=ES_NA_lines
#ES_data_wNA=cbind(ES_data,ES_NA_data)[,cell_lines]
#ES_data_complete=cbind(ES_var_names,ES_data_wNA)


### manually curated mutation data ###
MUTdatafile="/Users/adaemen/Dropbox/drug_predictors/mutations/COSMICmutationData_LSSVM_extended_v2.csv"
MUT_data_import=read.csv(MUTdatafile, row.names=1)
MUT_data=t(MUT_data_import)
MUT_var_names=cbind(rep("MM",length(rownames(MUT_data))),rownames(MUT_data))
colnames(MUT_var_names)=c("DataType","ID")

#Create NA columns for cell lines not represented
MUT_NA_lines=cell_lines[which(!cell_lines %in% colnames(MUT_data))]
MUT_NA_data=matrix(NA,nrow=length(rownames(MUT_data)), ncol=length(MUT_NA_lines))
colnames(MUT_NA_data)=MUT_NA_lines
MUT_data_wNA=cbind(MUT_data,MUT_NA_data)[,cell_lines]
MUT_data_complete=cbind(MUT_var_names,MUT_data_wNA)

### Amplification (GISTIC) data ###
#GISTdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/CNV/GISTIC_ROI_10JunV2.clean.txt"
#GIST_data_import=read.table(GISTdatafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)
#GIST_data=t(GIST_data_import)
#GIST_var_names=cbind(rep("CNV",length(rownames(GIST_data))),rownames(GIST_data))
#colnames(GIST_var_names)=c("DataType","ID")

#Create NA columns for cell lines not represented
#GIST_NA_lines=cell_lines[which(!cell_lines %in% colnames(GIST_data))]
#GIST_NA_data=matrix(NA,nrow=length(rownames(GIST_data)), ncol=length(GIST_NA_lines))
#colnames(GIST_NA_data)=GIST_NA_lines
#GIST_data_wNA=cbind(GIST_data,GIST_NA_data)[,cell_lines]
#GIST_data_complete=cbind(GIST_var_names,GIST_data_wNA)


### Amplification (SNP6) data ###
SNPdatafile="/Users/adaemen/Dropbox/drug_predictors/CNV/filtered_v2/SNP6_genelevel_stringent_std0.7.csv"
SNP_data_import=read.csv(SNPdatafile)
SNP_data=SNP_data_import[,6:length(colnames(SNP_data_import))]
SNP_feat_data=as.vector(SNP_data_import[,5])

SNP_var_names=cbind(rep("SNP",length(rownames(SNP_data))),SNP_feat_data)
colnames(SNP_var_names)=c("DataType","ID")

#Fix misnamed libraries
colnames(SNP_data)[which(colnames(SNP_data)=="X185A1")]="184A1"
colnames(SNP_data)[which(colnames(SNP_data)=="X184B5")]="184B5"
colnames(SNP_data)[which(colnames(SNP_data)=="X600MPE")]="600MPE"
colnames(SNP_data)[which(colnames(SNP_data)=="X21MT1")]="21MT1"
colnames(SNP_data)[which(colnames(SNP_data)=="X21MT2")]="21MT2"
colnames(SNP_data)[which(colnames(SNP_data)=="X21NT")]="21NT"
colnames(SNP_data)[which(colnames(SNP_data)=="X21PT")]="21PT"

#Create NA columns for cell lines not represented
SNP_NA_lines=cell_lines[which(!cell_lines %in% colnames(SNP_data))]
SNP_NA_data=matrix(NA,nrow=length(rownames(SNP_data)), ncol=length(SNP_NA_lines))
colnames(SNP_NA_data)=SNP_NA_lines
SNP_data_wNA=cbind(SNP_data,SNP_NA_data)[,cell_lines]
SNP_data_complete=cbind(SNP_var_names,SNP_data_wNA)


### RPPA data ###
RPPAdatafile="/Users/adaemen/Dropbox/drug_predictors/RPPA/filtered_v2/RPPA_validated.txt"
RPPA_data_import=read.table(RPPAdatafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)
RPPA_data=t(RPPA_data_import)
RPPA_var_names=cbind(rep("RPPA",length(rownames(RPPA_data))),rownames(RPPA_data))
colnames(RPPA_var_names)=c("DataType","ID")

#Fix misnamed libraries
colnames(RPPA_data)[which(colnames(RPPA_data)=="SUM102")]="SUM102PT"

#Create NA columns for cell lines not represented
RPPA_NA_lines=cell_lines[which(!cell_lines %in% colnames(RPPA_data))]
RPPA_NA_data=matrix(NA,nrow=length(rownames(RPPA_data)), ncol=length(RPPA_NA_lines))
colnames(RPPA_NA_data)=RPPA_NA_lines
RPPA_data_wNA=cbind(RPPA_data,RPPA_NA_data)[,cell_lines]
RPPA_data_complete=cbind(RPPA_var_names,RPPA_data_wNA)


### Western data ###
#Westdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/Western/JWGray_BCCL_Western_v2.txt"
#West_data_import=read.table(Westdatafile, header = TRUE, na.strings = "NA", sep="\t")
#West_feat_data=as.vector(West_data_import[,1])
#West_var_names=cbind(rep("WB",length(West_feat_data)),West_feat_data)
#colnames(West_var_names)=c("DataType","ID")
#West_data=West_data_import[,2:length(colnames(West_data_import))]
#West_data=log2(West_data+1) #Log2 data

#Fix misnamed libraries
#colnames(West_data)[which(colnames(West_data)=="X184A1")]="184A1"
#colnames(West_data)[which(colnames(West_data)=="X184B5")]="184B5"
#colnames(West_data)[which(colnames(West_data)=="X600MPE")]="600MPE"

#Create NA columns for cell lines not represented
#West_NA_lines=cell_lines[which(!cell_lines %in% colnames(West_data))]
#West_NA_data=matrix(NA,nrow=length(rownames(West_data)), ncol=length(West_NA_lines))
#colnames(West_NA_data)=West_NA_lines
#West_data_wNA=cbind(West_data,West_NA_data)[,cell_lines]
#West_data_complete=cbind(West_var_names,West_data_wNA)


### U133A data ###
U133Adatafile="/Users/adaemen/Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyRMA_genelevel_maxvar_stringent.csv"
U133A_data_import=read.csv(U133Adatafile, row.names=1)
U133A_var_names=cbind(rep("U133A",length(rownames(U133A_data_import))),rownames(U133A_data_import))
colnames(U133A_var_names)=c("DataType","ID")
U133A_data=U133A_data_import

#Fix misnamed libraries
colnames(U133A_data)[which(colnames(U133A_data)=="X600MPE")]="600MPE"

#Create NA columns for cell lines not represented
U133A_NA_lines=cell_lines[which(!cell_lines %in% colnames(U133A_data))]
U133A_NA_data=matrix(NA,nrow=length(rownames(U133A_data)), ncol=length(U133A_NA_lines))
colnames(U133A_NA_data)=U133A_NA_lines
U133A_data_wNA=cbind(U133A_data,U133A_NA_data)[,cell_lines]
U133A_data_complete=cbind(U133A_var_names,U133A_data_wNA)


### Methylation data ###
Methdatafile="/Users/adaemen/Dropbox/drug_predictors/methylation/filtered_v2/Methylation_stringent.csv"
MethAnnofile="/Users/adaemen/Dropbox/drug_predictors/methylation/filtered_v2/Methylation_annotation_stringent.csv"
Meth_data_import=read.csv(Methdatafile, row.names=1)
Meth_anno_data=read.csv(MethAnnofile)
#Create single combined ID: probe__Symbol
comb_ID=paste(Meth_anno_data[,1],Meth_anno_data[,4],sep="__")
Meth_var_names=cbind(rep("Meth",length(comb_ID)),comb_ID)
colnames(Meth_var_names)=c("DataType","ID")
Meth_data=Meth_data_import

#Fix misnamed libraries
colnames(Meth_data)[which(colnames(Meth_data)=="X600MPE")]="600MPE"

#Create NA columns for cell lines not represented
Meth_NA_lines=cell_lines[which(!cell_lines %in% colnames(Meth_data))]
Meth_NA_data=matrix(NA,nrow=length(rownames(Meth_data)), ncol=length(Meth_NA_lines))
colnames(Meth_NA_data)=Meth_NA_lines
Meth_data_wNA=cbind(Meth_data,Meth_NA_data)[,cell_lines]
Meth_data_complete=cbind(Meth_var_names,Meth_data_wNA)


#Create new dataframe/matrix for complete combined data
all_data=rbind(RNA_data_complete, EA_data_complete, U133A_data_complete, Meth_data_complete, MUT_data_complete, RPPA_data_complete, SNP_data_complete)
write.table(all_data, file=outfile, row.names=FALSE, quote = FALSE, sep="\t")
write.table(all_data, file=outfile2, row.names=FALSE, quote = FALSE, sep="\t", na="")



