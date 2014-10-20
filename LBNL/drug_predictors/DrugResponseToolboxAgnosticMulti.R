#!/usr/bin/env Rscript

# Copyright 2011 Anneleen Daemen, Obi Griffith, Laura Heiser

# DrugResponseToolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# DrugResponseToolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with DrugResponseToolbox. If not, see <http://www.gnu.org/licenses/>.

###Load raw data at individual patient-level (any combination of following data formats):
# Argument 1 = patient ID (string)
# Argument 2 = 'NA' or 1) Gene Expression data
# Argument 3 = 'NA' or 2) Methylation data. The values represent the proportion of methylated DNA at each single CpG locus.
# Argument 4 = 'NA' or 3) Copy Number data

### Example run on command line:
### ./DrugResponseToolboxAgnostic.R "FullCombination_test" "Expression_test_tumor.txt" "Meth_test_tumor.txt" "CNV_test_tumor.txt"

###Libraries
library(preprocessCore)
library(randomForest)

#args=c('FullCombination_test','Expression_test_tumor.txt','Meth_test_tumor.txt','CNV_test_tumor.txt')
#args=(commandArgs(TRUE))
#args.all=commandArgs(trailingOnly = FALSE)

###Determine path and name of script being run
#R_script=sub("--file=", "", args.all[grep("--file=", args.all)])
#R_script_name=strsplit(R_script,"/")[[1]][length(strsplit(R_script,"/")[[1]])] 
#script_path=paste(getwd(),R_script_name,sep="/")
#script_path=gsub("_",x=script_path,replacement="\\\\_") #Escapes "_" in path names

###Set work directory
#path='/Users/adaemen/Dropbox/drug_predictors/Rscripts/Rtoolbox/'
#path='C:/Users/Obi/Documents/My Dropbox/drug_predictors/Rscripts/Rtoolbox/'
#path='/csb/home/obig/Projects/drug_predictors/Rtoolbox/'
path='/Users/ogriffit/Dropbox/drug_predictors/Rscripts/Rtoolbox/'
setwd(path)
#pathPatientFiles='C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/TumorData/'
#pathPatientFiles='/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/TumorData/'
pathPatientFiles='/Users/ogriffit/Dropbox/drug_predictors/validation/RtoolboxResults/Breast/TumorData/'

#Load patient data file 
#Prevent NA from being intepretted as NA, since script expects this as string
#patientdata=read.table("C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/TCGA_agilent_meth_SNP6_308overlap_samples.txt",sep="\t",header=TRUE, as.is=c(1:4), na.strings = "") 
#patientdata=read.table("/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/CCLE_affymetrix_SNP6_samples.txt",sep="\t",header=TRUE, as.is=c(1:4), na.strings = "") 
patientdata=read.table("/Users/ogriffit/Dropbox/drug_predictors/validation/RtoolboxResults/Breast/BCCL_RNAseq_SNP6_samples.txt",sep="\t",header=TRUE, as.is=c(1:4), na.strings = "") 

###Results files (note tex will be converted to pdf with pdflatex
#Resultspath="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/Results/"
#Resultspath="/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/Results/"
#Resultspath="/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/Results_Fixed/"
Resultspath="/Users/ogriffit/Dropbox/drug_predictors/validation/RtoolboxResults/Breast/Results/"
#txt_output="TCGA_results_ExpMeth.txt" #For multi-patient run, summarize in single results file
#txt_output="TCGA_results_ExpMethCNV.txt" #For multi-patient run, summarize in single results file
#txt_output="CCLE_results_ExpCNV.txt" #For multi-patient run, summarize in single results file
txt_output="validation_results_ExpCNV.txt" #For multi-patient run, summarize in single results file

#Set which data types to use (to further restrict beyond what is available, e.g., to force use of expression/SNP only files)
#datatypes=c("Exp","Meth","CNV") #Use all data, where available
datatypes=c("Exp","CNV") #Use only expression and CNV, ignore Meth if present
#datatypes=c("Exp","Meth") #Use only expression and Meth, ignore CNV if present - testing purposes only

#Cutoff/Threshold values
WPMVcutoff=80 #Require at least X percent (weighted) of model variables to be available in input data
ModelAUCcutoff=0.6 #Require model training AUC of at least 0.6

for (i in 1:length(rownames(patientdata))){
#for (i in 1:2){ #testing purposes
	args=patientdata[i,]
	PatientID=args[1] #'FullCombination_test'
	Expfile=args[2] #'ExpressionData.txt'
	Methfile=args[3] #'MethData.txt'
	CNVfile=args[4] #'CNVData.txt'

	#Set file variables to NA if data type not desired
	if (!"Exp" %in% datatypes){Expfile='NA'}
	if (!"Meth" %in% datatypes){Methfile='NA'}
	if (!"CNV" %in% datatypes){CNVfile='NA'}

	if (Methfile=='NA'&CNVfile=='NA'&Expfile=='NA') {next} else 
	if (Methfile=='NA'&CNVfile=='NA') {DataCombi=1; DataCombiName="Expression"} else 
	if (Expfile=='NA'&CNVfile=='NA') {DataCombi=2; DataCombiName="Methylation"} else 
	if (Expfile=='NA'&Methfile=='NA') {DataCombi=3; DataCombiName="CNV"} else 
	if (CNVfile=='NA') {DataCombi=4; DataCombiName="Expression+Meth"} else 
	if (Methfile=='NA') {DataCombi=5; DataCombiName="Expression+CNV"} else 
	if (Expfile=='NA') {DataCombi=6; DataCombiName="Meth+CNV"} else 
	{DataCombi=7; DataCombiName="Expression+Meth+CNV"} 

	#Print current patient to screen to follow progress
	print(PatientID)

	###Load patient data and preprocess
	# Expression Data
	if (DataCombi==1|DataCombi==4|DataCombi==5|DataCombi==7) {
		# Get patient expression data
		ExpData=read.table(paste(pathPatientFiles,Expfile,sep=""),sep="\t",header=TRUE,row.names=1)
	
		# Normalize/scale patient data to cell line data (U133A, exon-array, RNAseq)
		#First load cell line data
		CellLineEAData=read.csv(paste(path,'ExtraFiles/breastExon_genelevel_stringent.csv',sep=""),sep=",",header=TRUE, row.names=1)
		CellLineU133AData=read.csv(paste(path,'ExtraFiles/Neve_AffyRMA_genelevel_maxvar_stringent.csv',sep=""),sep=",",header=TRUE, row.names=1)
		CellLineRSData=read.table(paste(path,'ExtraFiles/breastRNAseq_genelevel_stringent.txt',sep=""),sep="\t",header=TRUE, row.names=3)
		#Create mapping for RS data from gene name to long-format variable name, NOTE: this only works because in RNAseq stringent gene file, there is no ambiguity at gene symbol level
		RSmapping=as.data.frame(paste(CellLineRSData[,"FID"],rownames(CellLineRSData),rownames(CellLineRSData),sep="__"))
		rownames(RSmapping)=rownames(CellLineRSData)
		#Strip extra columns from RS data
		CellLineRSData=CellLineRSData[,4:length(colnames(CellLineRSData))]

		#Limit input data to only genes actually in cell line data, if gene missing from input data, set to NA
		ExpDataEA=ExpData[rownames(CellLineEAData),,drop=FALSE]
		rownames(ExpDataEA)=rownames(CellLineEAData)
		ExpDataU133A=ExpData[rownames(CellLineU133AData),,drop=FALSE]
		rownames(ExpDataU133A)=rownames(CellLineU133AData)
		ExpDataRS=ExpData[rownames(CellLineRSData),,drop=FALSE]
		rownames(ExpDataRS)=as.vector(RSmapping[,1]) #replace RS gene names with long format

		#Determine target distribution based on cell line data
		EA.target=normalize.quantiles.determine.target(as.matrix(CellLineEAData))
		RS.target=normalize.quantiles.determine.target(as.matrix(CellLineRSData))
		U133A.target=normalize.quantiles.determine.target(as.matrix(CellLineU133AData))

		#Set patient tumor data to target
		data.EA=normalize.quantiles.use.target(as.matrix(ExpDataEA),EA.target)
		data.RS=normalize.quantiles.use.target(as.matrix(ExpDataRS),RS.target)
		data.U133A=normalize.quantiles.use.target(as.matrix(ExpDataU133A),U133A.target)

		#NOTE - FOR RUNNING ON BCCL DATA, ABOVE STEP NOT NECESSARY - JUST USE EXISTING DATA
		data.EA=ExpDataEA
		data.RS=ExpDataRS
		data.U133A=ExpDataU133A
	
		#Set variable names to match those in models, and add patient ID as column heading
		data.EAG=data.EA #Necessary because different convention used for LSSVM and RF for naming of gene level EA variables
		data.RS2=data.RS #Necessary because different convention used for LSSVM and RF for naming of gene level variables (RS_ vs RS__)
		rownames(data.EA)=paste('EA__',rownames(ExpDataEA),sep="")
		rownames(data.EAG)=paste('EA_G__',rownames(ExpDataEA),sep="")
		rownames(data.RS)=paste('RS__',rownames(ExpDataRS),sep="")
		rownames(data.RS2)=paste('RS_',rownames(ExpDataRS),sep="")
		rownames(data.U133A)=paste('U133A__',rownames(ExpDataU133A),sep="")
		colnames(data.EA)=PatientID
		colnames(data.EAG)=PatientID
		colnames(data.RS)=PatientID
		colnames(data.RS2)=PatientID
		colnames(data.U133A)=PatientID
		setwd(path)
	}

	# Methylation - check values are within [0,1]
	if (DataCombi==2|DataCombi==4|DataCombi==6|DataCombi==7) {
		# Get patient expression data
		MethData=read.table(paste(pathPatientFiles,Methfile,sep=""),sep="\t",header=TRUE,row.names=1)

		# Normalize/scale patient data to cell line data
		#First load cell line data
		CellLineMethData=read.csv(paste(path,'ExtraFiles/Methylation_stringent.csv',sep=""),sep=",",header=TRUE, row.names=1)

		# Name conversion of cg-probe ids to genes
		#Limit input data to only genes actually in cell line data, if gene missing from input data, set to NA
		GeneProbeMapping.meth=read.csv(paste(path,'ExtraFiles/Methylation_annotation_stringent.csv',sep=""),sep=",",header=TRUE)
		MethData=as.data.frame(MethData[as.vector(GeneProbeMapping.meth[,1]),1])
		rownames(MethData)=GeneProbeMapping.meth[,1]

		#Determine target distribution based on cell line data
		Meth.target=normalize.quantiles.determine.target(as.matrix(CellLineMethData))

		#Set patient tumor data to target
		data.Meth=normalize.quantiles.use.target(as.matrix(MethData),Meth.target)

		#NOTE - FOR RUNNING ON BCCL DATA, ABOVE STEP NOT NECESSARY - JUST USE EXISTING DATA
		data.Meth=MethData

		#Set variable names to match those in models, and add patient ID as column heading
		rownames(data.Meth)=paste('Meth__',rownames(MethData),'__',as.vector(GeneProbeMapping.meth[,4]),sep="")
		colnames(data.Meth)=PatientID
		setwd(path)
	}

	# CNV data at gene level
	if (DataCombi==3|DataCombi==5|DataCombi==6|DataCombi==7) {
		# Get patient CNV data
		#CNVData=read.table(paste(pathPatientFiles,CNVfile,sep=""),sep="\t",header=TRUE, row.names=1)
		CNVData_raw=read.table(paste(pathPatientFiles,CNVfile,sep=""),sep="\t",header=TRUE)

		# Normalize/scale patient data to cell line data
		#First load cell line data
		CellLineCNVData=read.csv(paste(path,'ExtraFiles/SNP6/SNP6_genelevel_stringent_std0.7.csv',sep=""),sep=",",header=TRUE,row.names=5)
		CellLineCNVData=CellLineCNVData[,5:length(colnames(CellLineCNVData))]

		#Limit input data to only genes actually in cell line data, if gene missing from input data, set to NA
		overlap=which(CNVData_raw[,1] %in% rownames(CellLineCNVData))
		CNVData=as.data.frame(CNVData_raw[overlap,2,drop=FALSE])
		rownames(CNVData)=CNVData_raw[overlap,1]
		CNVData=CNVData[rownames(CellLineCNVData),,drop=FALSE]
		rownames(CNVData)=rownames(CellLineCNVData)

		#Determine target distribution based on cell line data
		CNV.target=normalize.quantiles.determine.target(as.matrix(CellLineCNVData))

		#Set patient tumor data to target
		data.SNP=normalize.quantiles.use.target(as.matrix(CNVData),CNV.target)

		#NOTE - FOR RUNNING ON BCCL DATA, ABOVE STEP NOT NECESSARY - JUST USE EXISTING DATA
		data.SNP=CNVData

		#Set variable names to match those in models, and add patient ID as column heading
		rownames(data.SNP)=paste('SNP__',rownames(CNVData),sep="")
		colnames(data.SNP)=PatientID
		setwd(path)
	}

	setwd(Resultspath)

	### Testing of best predictor per drug compound on tumor sample (best predictor defined upon the inputed data combination) for set of 91 drug compounds
	drugs_interest=read.table(paste(path,'ExtraFiles/DrugCompounds.txt',sep=""),header=TRUE,sep="\t")
	predictor_perf = data.frame(cbind(Sample=NA, InputData=NA, drugs_interest, Model=NA, ModelAUC=NA, DataCombination=NA, NbModelVars=NA, PercModelVars=NA, WtPercModelVars=NA, Score=NA, Probability=NA, Sensitivity=NA), stringsAsFactors=FALSE, row.names=drugs_interest[[1]])

	### Mapping table with best predictor per drug compound
	### Drug x Data Combi matrix with best predictor given a data combi (might be based on only 1 data set, despite of available other data types)
	#predictor_mapping=read.table(paste(path,'ExtraFiles/PredictorsMapping_PlatformAgnosticToolbox_Expression.txt',sep=""),sep="\t",header=TRUE,row.names=1)
	#predictor_mapping_AUC=read.table(paste(path,'ExtraFiles/PredictorsMapping_PlatformAgnosticToolbox_Expression_AUC.txt',sep=""),sep="\t",header=TRUE,row.names=1)

	predictor_mapping=read.table(paste(path,'ExtraFiles/PredictorsMapping_PlatformAgnosticToolbox_RNAseq.txt',sep=""),sep="\t",header=TRUE,row.names=1)
	predictor_mapping_AUC=read.table(paste(path,'ExtraFiles/PredictorsMapping_PlatformAgnosticToolbox_RNAseq_AUC.txt',sep=""),sep="\t",header=TRUE,row.names=1)


	for (drug_compound in drugs_interest[[1]]) {
		#drug_compound=as.vector(drugs_interest[[1]][40]) #for testing purposes
		drug_compound=as.vector(drug_compound)
		predictor_perf[drug_compound,"Sample"]=PatientID
		predictor_perf[drug_compound,"InputData"]=DataCombiName
		model=as.character(predictor_mapping[drug_compound,DataCombi])
		modelAUC=predictor_mapping_AUC[drug_compound,DataCombi]
		print(paste("processing ",drug_compound," with ", model))
		if (regexpr('_U133A_Meth_SNP_',model)[1]>0) {
			DataCombiModel="U133A-Meth-SNP"
			CombiData_LSSVM=rbind(data.U133A,data.Meth,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianU133A=median(data.U133A[!is.na(data.U133A)])
			data.U133A[is.na(data.U133A)]=medianU133A
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.U133A,data.Meth,data.SNP)
		} else if (regexpr('_RS_Meth_SNP_',model)[1]>0) {
			DataCombiModel="RS-Meth-SNP"
			CombiData_LSSVM=rbind(data.RS2,data.Meth,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.RS2)[1]),rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianRS=median(data.RS[!is.na(data.RS)])
			data.RS[is.na(data.RS)]=medianRS
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.RS,data.Meth,data.SNP)
		} else if (regexpr('_EA_Meth_SNP_',model)[1]>0) {
			DataCombiModel="EA-Meth-SNP"
			CombiData_LSSVM=rbind(data.EAG,data.Meth,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.EAG)[1]),rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianEA=median(data.EA[!is.na(data.EA)])
			data.EA[is.na(data.EA)]=medianEA
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.EA,data.Meth,data.SNP)
		} else if (regexpr('_U133A_Meth_',model)[1]>0) {
			DataCombiModel="U133A-Meth"
			CombiData_LSSVM=rbind(data.U133A,data.Meth)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(1,dim(data.Meth)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianU133A=median(data.U133A[!is.na(data.U133A)])
			data.U133A[is.na(data.U133A)]=medianU133A
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			CombiData_RF=rbind(data.U133A,data.Meth)
		} else if (regexpr('_U133A_SNP_',model)[1]>0) {
			DataCombiModel="U133A-SNP"
			CombiData_LSSVM=rbind(data.U133A,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianU133A=median(data.U133A[!is.na(data.U133A)])
			data.U133A[is.na(data.U133A)]=medianU133A
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.U133A,data.SNP)
		} else if (regexpr('_RS_Meth_',model)[1]>0) {
			DataCombiModel="RS-Meth"
			CombiData_LSSVM=rbind(data.RS2,data.Meth)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.RS2)[1]),rep(1,dim(data.Meth)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianRS=median(data.RS[!is.na(data.RS)])
			data.RS[is.na(data.RS)]=medianRS
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			CombiData_RF=rbind(data.RS,data.Meth)
		} else if (regexpr('_RS_SNP_',model)[1]>0) {
			DataCombiModel="RS-SNP"
			CombiData_LSSVM=rbind(data.RS2,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.RS2)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianRS=median(data.RS[!is.na(data.RS)])
			data.RS[is.na(data.RS)]=medianRS
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.RS,data.SNP)
		} else if (regexpr('_EA_Meth_',model)[1]>0) {
			DataCombiModel="EA-Meth"
			CombiData_LSSVM=rbind(data.EAG,data.Meth)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.EAG)[1]),rep(1,dim(data.Meth)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianEA=median(data.EA[!is.na(data.EA)])
			data.EA[is.na(data.EA)]=medianEA
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			CombiData_RF=rbind(data.EA,data.Meth)
		} else if (regexpr('_EA_SNP_',model)[1]>0) {
			DataCombiModel="EA-SNP"
			CombiData_LSSVM=rbind(data.EAG,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.EAG)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianEA=median(data.EA[!is.na(data.EA)])
			data.EA[is.na(data.EA)]=medianEA
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.EA,data.SNP)
		} else if (regexpr('_Meth_SNP_',model)[1]>0) {
			DataCombiModel="Meth-SNP"
			CombiData_LSSVM=rbind(data.Meth,data.SNP)
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(1,dim(data.Meth)[1]),rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=rbind(data.Meth,data.SNP)
		} else if (regexpr('_U133A_',model)[1]>0) {
			DataCombiModel="U133A"
			CombiData_LSSVM=data.U133A
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.U133A)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianU133A=median(data.U133A[!is.na(data.U133A)])
			data.U133A[is.na(data.U133A)]=medianU133A
			CombiData_RF=data.U133A
		} else if (regexpr('_RS_',model)[1]>0) {
			DataCombiModel="RS"
			CombiData_LSSVM=data.RS2
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.RS2)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianRS=median(data.RS[!is.na(data.RS)])
			data.RS[is.na(data.RS)]=medianRS
			CombiData_RF=data.RS
		} else if (regexpr('_EA_',model)[1]>0) {
			DataCombiModel="EA"
			CombiData_LSSVM=data.EAG
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(11,dim(data.EAG)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianEA=median(data.EA[!is.na(data.EA)])
			data.EA[is.na(data.EA)]=medianEA
			CombiData_RF=data.EA
		} else if (regexpr('_Meth_',model)[1]>0) {
			DataCombiModel="Meth"
			CombiData_LSSVM=data.Meth
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(1,dim(data.Meth)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianMeth=median(data.Meth[!is.na(data.Meth)])
			data.Meth[is.na(data.Meth)]=medianMeth
			CombiData_RF=data.Meth
		} else if (regexpr('_SNP_',model)[1]>0) {
			DataCombiModel="SNP"
			CombiData_LSSVM=data.SNP
			rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
			rangeClKernel[1,]=c(rep(10,dim(data.SNP)[1]))
			colnames(rangeClKernel)=rownames(CombiData_LSSVM)
			medianSNP=median(data.SNP[!is.na(data.SNP)])
			data.SNP[is.na(data.SNP)]=medianSNP
			CombiData_RF=data.SNP
		}
	
		if (regexpr('LSSVM',model)[1]>0) {
			ModelDrug=read.table(paste(path,'ModelsAgnostic/',model,'.txt',sep=""),as.is=1)
			TrainDrug=read.table(paste(path,'ModelsAgnostic/Trainingdata_',model,'.txt',sep=""),header=TRUE,sep="\t",row.names=1)

			# Compare input data variables to model variables and determine if enough present for prediction
			model_variables=rownames(TrainDrug)[2:length(rownames(TrainDrug))]
			#First calculate simple percent of model variables found in the input data
			PercModelVars=(length(which(rownames(CombiData_LSSVM)[which(!is.na(CombiData_LSSVM))] %in% model_variables))/length(model_variables))*100
			#Next calculate a weighted percent where higher ranking variables count more towards percent
			model_variable_ranks=seq(from=length(model_variables), to=1) #Best variable = largest rank value. LSSVM model variables already ordered by rank
			input_variable_ranks=model_variable_ranks[which(model_variables %in% rownames(CombiData_LSSVM)[which(!is.na(CombiData_LSSVM))])]
			WtPercModelVars=(sum(input_variable_ranks)/sum(model_variable_ranks))*100

			Asigm=ModelDrug[[1]][1]
			Bsigm=ModelDrug[[1]][2]
			bcoeff=ModelDrug[[1]][3]
			alpha=ModelDrug[[1]][-(1:3)]
			target=as.matrix(TrainDrug[1,])
			traindata=as.matrix(TrainDrug[-1,])
			FeatureSubset=rownames(TrainDrug[-1,])
			CombiSubData=as.data.frame(CombiData_LSSVM[FeatureSubset,])
			rownames(CombiSubData)=FeatureSubset
			rangeClKernelSub=rangeClKernel[1,FeatureSubset]
		
			### Reduction to model features that are present in the input data file(s)
			removal=which(is.na(CombiSubData[,1]))
			if (length(removal)>0) {	
				traindata=traindata[-removal,]
				FeatureSubset=FeatureSubset[-removal]
				CombiSubData=as.data.frame(CombiSubData[-removal,])
				rownames(CombiSubData)=FeatureSubset
				rangeClKernel=rangeClKernelSub[FeatureSubset]
			}
		
			Kfeature=0
			clinicalKernel <- function(i) {
				clinicalKernelIntern <- function(j) {
					Kfeature[j]=(rangeClKernel[i]-abs(traindata[i,j]-CombiSubData[i,]))
				}
				Kfeature=sapply(c(1:dim(traindata)[2]),clinicalKernelIntern)
				Kfeature=Kfeature/rangeClKernel[i]
			}
			KtestPerFeature=sapply(c(1:length(FeatureSubset)),clinicalKernel)
			Ktest=rowSums(KtestPerFeature)/length(FeatureSubset)
			#Ktest=t(traindata) %*% as.matrix(CombiSubData)
			latent=((alpha*target) %*% Ktest) + bcoeff
			predictor_perf[drug_compound,"Score"]=format(latent,digits=5)
			prob=1/(1+exp(Asigm*latent+Bsigm))
			predictor_perf[drug_compound,"Probability"]=format(prob,digits=5)
			predictor_perf[drug_compound,"Sensitivity"]=(prob>0.5)
			predictor_perf[drug_compound,"Model"]="LSSVM"
			predictor_perf[drug_compound,"DataCombination"]=DataCombiModel
			predictor_perf[drug_compound,"PercModelVars"]=PercModelVars
			predictor_perf[drug_compound,"WtPercModelVars"]=WtPercModelVars
			predictor_perf[drug_compound,"NbModelVars"]=length(model_variables)
			predictor_perf[drug_compound,"ModelAUC"]=modelAUC
		} else if (regexpr('RF',model)[1]>0) {
			RFdata=t(CombiData_RF)
			load(file=paste(path,'ModelsAgnostic/',model,'.Rdata',sep=""))
	
			# Compare input data variables to model variables and determine if enough present for prediction
			# use CombiData_LSSVM to determing NA values because they were replaced with imputed values for CombiData_RF
			model_variables=names(rf_model$importance[,1])
			variable_importances=rf_model$importance[,4]
			#First calculate simple percent of model variables found in the input data
			PercModelVars=(length(which(rownames(CombiData_RF)[which(!is.na(CombiData_LSSVM))] %in% model_variables))/length(model_variables))*100
			#Next calculate a weighted percent where higher ranking variables count more towards percent
			model_variable_ranks=order(variable_importances) #Best variable (largest var imp) = largest rank value
			input_variable_ranks=model_variable_ranks[which(model_variables %in% rownames(CombiData_RF)[which(!is.na(CombiData_LSSVM))])]
			WtPercModelVars=(sum(input_variable_ranks)/sum(model_variable_ranks))*100

			RF_predictions_response=predict(rf_model, RFdata, type="response")
			RF_predictions_prob=predict(rf_model, RFdata, type="prob")
			RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
			predictor_perf[drug_compound,"Score"]=RF_predictions_vote[1,"sensitive"]
			predictor_perf[drug_compound,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
			predictor_perf[drug_compound,"Sensitivity"]=RF_predictions_response=="sensitive"
			predictor_perf[drug_compound,"Model"]="RF"
			predictor_perf[drug_compound,"DataCombination"]=DataCombiModel
			predictor_perf[drug_compound,"PercModelVars"]=PercModelVars
			predictor_perf[drug_compound,"WtPercModelVars"]=WtPercModelVars
			predictor_perf[drug_compound,"NbModelVars"]=length(model_variables)
			predictor_perf[drug_compound,"ModelAUC"]=modelAUC
		}
	}

	### Convert Sensitivity=TRUE/FALSE to +/-
	predictor_perf[which(predictor_perf[,"Sensitivity"]==TRUE),"Sensitivity"]="+"
	predictor_perf[which(predictor_perf[,"Sensitivity"]==FALSE),"Sensitivity"]="-"

	#Apply cutoffs to limit results to only those drugs with sufficient input data
	predictor_perf[predictor_perf[,"WtPercModelVars"]<WPMVcutoff,"Score"]=NA
	predictor_perf[predictor_perf[,"WtPercModelVars"]<WPMVcutoff,"Probability"]=NA
	predictor_perf[predictor_perf[,"WtPercModelVars"]<WPMVcutoff,"Sensitivity"]=NA
	predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Score"]=NA
	predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Probability"]=NA
	predictor_perf[predictor_perf[,"ModelAUC"]<ModelAUCcutoff,"Sensitivity"]=NA

	if (i==1){
		predictor_perf_all=predictor_perf
	}else{
		predictor_perf_all=rbind(predictor_perf_all,predictor_perf)
	}
}
### Save results to txt file
write.table(predictor_perf_all, file=txt_output, sep="\t", row.names=FALSE)

