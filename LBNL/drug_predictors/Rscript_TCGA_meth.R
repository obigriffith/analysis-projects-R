library(Biobase)
library(affy)
library(randomForest)
library(DNAcopy) #BioConductor
library(CNTools) #BioConductor

path='/Users/adaemen/Dropbox/drug_predictors/Rscripts/Rtoolbox/'
setwd(path)

TCGAmeth_samples=read.table(paste(path,'TumorData/TCGAmeth_samples.txt',sep=""),sep="\t",header=FALSE)

DrugResponseFunction <- function(i) {
	PatientID=TCGAmeth_samples[[1]][i]
	Methfile=paste("Meth_TCGA_",TCGAmeth_samples[[1]][i],".txt",sep="")
	DataCombi=2
	
	###Set work directory
	path='/Users/adaemen/Dropbox/drug_predictors/Rscripts/Rtoolbox/'
	setwd(path)
	pathPatientFiles='TumorData/'
	
	###Results files (note tex will be converted to pdf with pdflatex
	txt_output=paste(PatientID,"_results.txt",sep="")
	
	###Load patient data and preprocess
	# Methylation - check values are within [0,1]
	MethData=read.table(paste(pathPatientFiles,Methfile,sep=""),sep="\t")
	rownames(MethData)=MethData[,1]
	
	# Name conversion of cg-probe ids to genes
	GeneProbeMapping.meth=read.csv(paste(path,'ExtraFiles/Methylation_annotation_stringent.csv',sep=""),sep=",",header=TRUE)
	data.Meth=as.data.frame(MethData[as.vector(GeneProbeMapping.meth[,1]),2])
	rownames(data.Meth)=GeneProbeMapping.meth[,1]
	rownames(data.Meth)=paste('Meth__',rownames(data.Meth),'__',as.vector(GeneProbeMapping.meth[,4]),sep="")
	colnames(data.Meth)=PatientID
	setwd(paste(path,'Results/',sep=""))

	### Testing of best predictor per drug compound on tumor sample (best predictor defined upon the inputed data combination) for set of 91 drug compounds
	drugs_interest=read.table(paste(path,'ExtraFiles/DrugCompounds.txt',sep=""),header=TRUE,sep="\t")
	predictor_perf = data.frame(cbind(drugs_interest, Model=NA, DataCombination=NA, Score=NA, Probability=NA, Sensitivity=NA), stringsAsFactors=FALSE)
	rownames(predictor_perf)=drugs_interest[[1]]

	### Mapping table with best predictor per drug compound
	### Drug x Data Combi matrix with best predictor given a data combi (might be based on only 1 data set, despite of available other data types)
	predictor_mapping=read.table(paste(path,'ExtraFiles/PredictorsMapping.txt',sep=""),sep="\t",header=TRUE,row.names=1)

	for (drug_compound in drugs_interest[[1]]) {
		drug_compound=as.vector(drug_compound)
		model=as.character(predictor_mapping[drug_compound,DataCombi])
		CombiData_LSSVM=data.Meth
		medianMeth=median(data.Meth[!is.na(data.Meth)])
		data.Meth[is.na(data.Meth)]=medianMeth
		CombiData_RF=rbind(data.Meth)
		DataCombiModel="Meth"
		rangeClKernel=mat.or.vec(1,dim(CombiData_LSSVM)[1])
		rangeClKernel[1,]=c(rep(1,dim(data.Meth)[1]))
		colnames(rangeClKernel)=rownames(CombiData_LSSVM)
	
	if (regexpr('LSSVM',model)[1]>0) {
		ModelDrug=read.table(paste(path,'Models/',model,'.txt',sep=""),as.is=1)
		TrainDrug=read.table(paste(path,'Models/Trainingdata_',model,'.txt',sep=""),header=TRUE,sep="\t",row.names=1)
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
	} else if (regexpr('RF',model)[1]>0) {
		RFdata=t(CombiData_RF)
		load(file=paste(path,'Models/',model,'.RData',sep=""))
		RF_predictions_response=predict(rf_model, RFdata, type="response")
		RF_predictions_prob=predict(rf_model, RFdata, type="prob")
		RF_predictions_vote=predict(rf_model, RFdata, type="vote", norm.votes=FALSE)
		predictor_perf[drug_compound,"Score"]=RF_predictions_vote[1,"sensitive"]
		predictor_perf[drug_compound,"Probability"]=format(RF_predictions_prob[1,"sensitive"],digits=5)
		predictor_perf[drug_compound,"Sensitivity"]=RF_predictions_response=="sensitive"
		predictor_perf[drug_compound,"Model"]="RF"
		predictor_perf[drug_compound,"DataCombination"]=DataCombiModel
	}
	}	

	### Save results to txt file
	write.table(predictor_perf, file=txt_output, sep="\t", row.names=FALSE)
}

MethResults=sapply(c(1:dim(TCGAmeth_samples)[1]),DrugResponseFunction)
