#library(kernlab)
### LS-SVM {kernlab}
#lssvm_model=lssvm(t(data_nonNA), target, scaled=FALSE, kernel="vanilladot", tau=0.01, reduced=FALSE, na.action=na.omit)

library(gplots)
library(ROCR)
library(matlab)

# machine precision
eps=2.2e-16

# nb of randomizations
nbRandomizations=500

#Specify working directories and predictor data file (4 different versions) - choose one
#RMA - Standard CDF
outdir="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/RMA_standardCDF/"
datafile="/Users/adaemen/Dropbox/drug_predictors/U133A/filtered/Neve_AffyRMA_genelevel_maxvar_stringent.csv"
outfile="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/RMA_standardCDF/LSSVM_results_summary_RMA_standardCDF.txt"

#GCRMA - standard CDF
#outdir="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/GCRMA_standardCDF/"
#datafile="/Users/adaemen/Dropbox/drug_predictors/U133A/filtered/Neve_AffyGCRMA_genelevel_maxvar_stringent.csv"
#outfile="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/GCRMA_standardCDF/LSSVM_results_summary_GCRMA_standardCDF.txt"

#RMA - custom CDF
#outdir="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/RMA_customCDF/"
#datafile="/Users/adaemen/Dropbox/drug_predictors/U133A/filtered/Neve_AffyRMAcustom_genelevel_stringent.csv"
#outfile="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/RMA_customCDF/LSSVM_results_summary_RMA_customCDF.txt"

#GCRMA custom CDF
#outdir="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/GCRMA_customCDF/"
#datafile="/Users/adaemen/Dropbox/drug_predictors/U133A/filtered/Neve_AffyGCRMAcustom_genelevel_stringent.csv"
#outfile="/Users/adaemen/Dropbox/drug_predictors/LSSVMpredictors/U133A/DrugsOfInterest/GCRMA_customCDF/LSSVM_results_summary_GCRMA_customCDF.txt"

#Set working directory and output files
setwd(outdir)

#Import combined predictor data
raw_data_import=read.csv(datafile, row.names=1)

#Break data into features info and expression values
raw_feat_data=rownames(raw_data_import)
raw_data=raw_data_import

#Fix column names
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"

#Import cell line data
cell_line_datafile="/Users/adaemen/Dropbox/drug_predictors/BCCL_data_list.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Filter down to set of cell lines with both drug data and U133A data (exclude non-malignant, non-BCCL)
core_cell_lines=c("BT20","HCC1143","HCC1187","HCC1500","HCC1569","HCC1937","HCC1954","HCC3153","HCC70","MDAMB468","SUM149PT","BT549","HCC38","HS578T","MDAMB157","MDAMB231","MDAMB436","SUM1315MO2","SUM159PT","600MPE","AU565","BT474","BT483","CAMA1","HCC1428","HCC202","HCC2185","LY2","MCF7","MDAMB134VI","MDAMB175VII","MDAMB361","MDAMB415","MDAMB453","SKBR3","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","ZR751","ZR7530","ZR75B")
data=raw_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#drug response data
drugdatafile="/Users/adaemen/Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_10Feb.csv"
raw_drugdata_import=read.csv(drugdatafile)
drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Fix row names
rownames(drug_data)[which(rownames(drug_data)=="Hs578T")]="HS578T"

#Retrieve data for only libraries in core cell line set
drug_data_filt=drug_data[core_cell_lines,]

#Transform data to -log values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#Load mean GI50 values from file
meanGI50file="/Users/adaemen/Dropbox/drug_predictors/drugdata/GI50meanThresholds.csv"
meanGI50_import=read.csv(meanGI50file, row.names=1)
drugnames=rownames(meanGI50_import)
meanGI50_data=meanGI50_import[,1]

#Fix problem drug name:
drugnames[drugnames=="Tykerb:IGF1R (1:1)"]="Tykerb(IGF1R)"

# Drugs of interest
drugs_interest=c("X5.FU","AKT1.2.inhibitor","Cisplatin","Docetaxel","SAHA..Vorinostat.","Vinorelbine","GSK2126458A","GSK_Tykerb")

#For each drug of interest, find predictors associated with response and build predictor for drug response
#First, create dataframe to hold results
predictor_perf = data.frame(cbind(drugs_interest, N=NA, meanAUCtest=NA, stdAUCtest=NA, meanAUCtrain=NA, stdAUCtrain=NA, meanGI50=NA), stringsAsFactors=FALSE)
rownames(predictor_perf)=drugs_interest

for (drug in drugs_interest) {
	i=which(drugs==drug)
	drugname=drugnames[i]
	print(paste("processing",drugname))
	
	#Divide cell lines into responders and non-responders
	drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
	names(drug_data_interest)=rownames(drug_data_filt_trans)
	cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
	drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
	
	#Retrieve mean from pre-calculated file
	mean_cutoff=meanGI50_data[i]
	
	drug_data_interest_NA=which(is.na(drug_data_interest))
	resistants=which(drug_data_interest<=mean_cutoff)
	sensitives=which(drug_data_interest>mean_cutoff)
	
	response_class=vector(length=length(drug_data_interest))
	response_class[drug_data_interest_NA]=NA
	response_class[sensitives]="sensitive"
	response_class[resistants]="resistant"
	response_class2=vector(length=length(drug_data_interest))
	response_class2[drug_data_interest_NA]=NA
	response_class2[sensitives]=1
	response_class2[resistants]=-1
	
	#Exclude libs where response_class=NA
	nonNA=which(!is.na(response_class))
	data_nonNA=data[,nonNA]
	response_class_nonNA=response_class[nonNA]
	response_class2_nonNA=response_class2[nonNA]
	cell_line_data_nonNA=cell_line_data[nonNA,]
	
	target=as.factor(response_class_nonNA)
	target2=response_class2_nonNA
	
	
	################################
	### Implementation of LS-SVM ###
	################################
	# Randomization, with split of samples in 2/3rd training, 1/3rd testing, stratified to outcome
	RandomizationVector=c(1:nbRandomizations)
	AUCrand_train=matrix(data=NA,nrow=1,ncol=length(RandomizationVector))
	AUCrand_test=matrix(data=NA,nrow=1,ncol=length(RandomizationVector))
	GammaRand=matrix(data=NA,nrow=1,ncol=length(RandomizationVector))
	
	for (rand in RandomizationVector) {
		indices_sens_rand=which(response_class_nonNA=="sensitive")
		indices_res_rand=which(response_class_nonNA=="resistant")
		trainindices_sens=sample(indices_sens_rand,round(length(indices_sens_rand)*2/3))
		trainindices_res=sample(indices_res_rand,round(length(indices_res_rand)*2/3))
		testindices_sens=setdiff(indices_sens_rand,trainindices_sens)
		testindices_res=setdiff(indices_res_rand,trainindices_res)
		trainindices=c(trainindices_sens,trainindices_res)
		testindices=c(testindices_sens,testindices_res)
		
		TrainRand=as.matrix(data_nonNA[,trainindices])
		TestRand=as.matrix(data_nonNA[,testindices])
		traintarget_rand=target2[trainindices]
		testtarget_rand=target2[testindices]
		response_class_rand=response_class_nonNA[trainindices]
		
		# 5-fold training cross-validation, with split stratified to outcome
		nrFolds=5
		indices_sens=which(response_class_rand=="sensitive")
		indices_res=which(response_class_rand=="resistant")
		num_sensitive=length(indices_sens)
		num_resistant=length(indices_res)
		Ind_sens=rep(0,num_sensitive)
		folds=ceiling(nrFolds*c(1:num_sensitive)/num_sensitive)
		Kperm=sample(c(1:nrFolds),nrFolds)
		Nperm=sample(c(1:num_sensitive),num_sensitive)
		Ind_sens[Nperm]=Kperm[folds]
		Ind_res=rep(0,num_resistant)
		folds=ceiling(nrFolds*c(1:num_resistant)/num_resistant)
		Kperm=sample(c(1:nrFolds),nrFolds)
		Nperm=sample(c(1:num_resistant),num_resistant)
		Ind_res[Nperm]=Kperm[folds]
		
		# Result vector with latent variables
		GammaVector=logspace(-4,6,n=40)
		counter=0
		AUCvector=matrix(data=NA,nrow=1,ncol=length(GammaVector))
		
		for (gamma in GammaVector) {
			counter=counter+1
			LatentVector=matrix(data=NA,nrow=1,ncol=length(traintarget_rand))
			for (fold in c(1:5)) {
				testIndPos=which(Ind_sens==fold)
				testIndNeg=which(Ind_res==fold)
				trainIndPos=which(Ind_sens!=fold)
				trainIndNeg=which(Ind_res!=fold)
				Test=c(indices_sens[testIndPos],indices_res[testIndNeg])
				Train=c(indices_sens[trainIndPos],indices_res[trainIndNeg])
				
				traindata=as.matrix(TrainRand[,Train])
				testdata=as.matrix(TrainRand[,Test])
				traintarget=traintarget_rand[Train]
				testtarget=traintarget_rand[Test]
				
				# Kernel matrix based on (normalized) linear kernel function
				Ktrain=t(traindata) %*% traindata
				#Dtrain=diag(1/sqrt(diag(Xtrain)+eps))
				#Ktrain=Dtrain %*% Xtrain %*% Dtrain
				
				Xtraintest=t(traindata) %*% testdata
				#Xtest=t(testdata) %*% testdata
				#Dtest = diag(1/sqrt(diag(Xtest)+eps))	
				#Ktest = Dtrain %*% Xtraintest %*% Dtest
				Ktest=Xtraintest
				
				# Training of LS-SVM model
				dummy=traintarget %*% t(traintarget)
				Omega=dummy*Ktrain
				A=matrix(data=NA,nrow=length(traintarget)+1,ncol=length(traintarget)+1)
				A[,1]=c(0,traintarget)
				A[1,]=c(0,t(traintarget))
				A[-1,-1]=Omega+diag(length(traintarget))/gamma;
				B=0
				B=c(B,rep(1,length(traintarget)))
				
				X=solve(A,B)
				
				b=X[1]
				alpha=X[2:length(X)]
				
				# Testing of LS-SVM model
				prediction_latent=(t(alpha*traintarget) %*% Ktest) + b
				prediction_label=sign(prediction_latent-eps)
				LatentVector[Test]=prediction_latent
			}
			
			# Calculation of 5-fold CV test performance for each gamma
			pred=prediction(as.vector(LatentVector),traintarget_rand)
			perf_AUC=performance(pred,"auc")
			AUC=perf_AUC@y.values[[1]]
			AUCvector[counter]=AUC
		}
		
		# Selection of optimal gamma, corresponding to highest AUC
		AUCmax=max(AUCvector)
		index_max=which(AUCvector==AUCmax)
		
		# Second round to further refine gamma
		if (index_max[1]==1) {
			index_lower=index_max[1]
			index_upper=index_max[1]+1
		} else if (index_max[1]==40) {
			index_lower=index_max[1]-1
			index_upper=index_max[1]
		} else {
			index_lower=index_max[1]-1
			index_upper=index_max[1]+1
		}
		
		GammaVectorRefine=logspace(log10(GammaVector[index_lower]),log10(GammaVector[index_upper]),n=40)
		counter=0
		AUCvectorRefine=matrix(data=NA,nrow=1,ncol=length(GammaVector))
		
		for (gamma in GammaVectorRefine) {
			counter=counter+1
			LatentVector=matrix(data=NA,nrow=1,ncol=length(traintarget_rand))
			for (fold in c(1:5)) {
				testIndPos=which(Ind_sens==fold)
				testIndNeg=which(Ind_res==fold)
				trainIndPos=which(Ind_sens!=fold)
				trainIndNeg=which(Ind_res!=fold)
				Test=c(indices_sens[testIndPos],indices_res[testIndNeg])
				Train=c(indices_sens[trainIndPos],indices_res[trainIndNeg])
				
				traindata=as.matrix(TrainRand[,Train])
				testdata=as.matrix(TrainRand[,Test])
				traintarget=traintarget_rand[Train]
				testtarget=traintarget_rand[Test]
				
				# Kernel matrix based on (normalized) linear kernel function
				Ktrain=t(traindata) %*% traindata
				
				Ktest=t(traindata) %*% testdata
				
				# Training of LS-SVM model
				dummy=traintarget %*% t(traintarget)
				Omega=dummy*Ktrain
				A=matrix(data=NA,nrow=length(traintarget)+1,ncol=length(traintarget)+1)
				A[,1]=c(0,traintarget)
				A[1,]=c(0,t(traintarget))
				A[-1,-1]=Omega+diag(length(traintarget))/gamma;
				B=0
				B=c(B,rep(1,length(traintarget)))
				
				X=solve(A,B)
				
				b=X[1]
				alpha=X[2:length(X)]
				
				# Testing of LS-SVM model
				prediction_latent=(t(alpha*traintarget) %*% Ktest) + b
				prediction_label=sign(prediction_latent-eps)
				LatentVector[Test]=prediction_latent
			}
			
			# Calculation of 5-fold CV test performance for each gamma
			pred=prediction(as.vector(LatentVector),traintarget_rand)
			perf_AUC=performance(pred,"auc")
			AUC=perf_AUC@y.values[[1]]
			AUCvectorRefine[counter]=AUC
		}
		
		# Selection of optimal gamma, corresponding to highest AUC
		AUCmax=max(AUCvectorRefine)
		index_max=which(AUCvectorRefine==AUCmax)
		GammaOpt=GammaVectorRefine[index_max[1]]
		
		# Retrain model on full training set with optimal gamma
		Ktrain=t(TrainRand) %*% TrainRand
		Ktest=t(TrainRand) %*% TestRand
		
		dummy=traintarget_rand %*% t(traintarget_rand)
		Omega=dummy*Ktrain
		A=matrix(data=NA,nrow=length(traintarget_rand)+1,ncol=length(traintarget_rand)+1)
		A[,1]=c(0,traintarget_rand)
		A[1,]=c(0,t(traintarget_rand))
		A[-1,-1]=Omega+diag(length(traintarget_rand))/GammaOpt;
		B=0
		B=c(B,rep(1,length(traintarget_rand)))
		X=solve(A,B)
		b=X[1]
		alpha=X[2:length(X)]
		
		prediction_latent=(t(alpha*traintarget_rand) %*% Ktest) + b
		prediction_label=sign(prediction_latent-eps)
		LatentVectorRand=prediction_latent
		
		pred=prediction(as.vector(LatentVectorRand),testtarget_rand)
		perf_AUC=performance(pred,"auc")
		AUCrand=perf_AUC@y.values[[1]]
		AUCrand_test[rand]=AUCrand
		AUCrand_train[rand]=AUCmax
		GammaRand[rand]=GammaOpt
	} #end of randomization loop
	
	mean(AUCrand_test)
	std(t(AUCrand_test))
	mean(AUCrand_train)
	std(t(AUCrand_train))
	
	#Add performance stats to dataframe
	predictor_perf[drug,"N"]=length(target2)
	predictor_perf[drug,"meanAUCtest"]=mean(AUCrand_test)
	predictor_perf[drug,"stdAUCtest"]=std(t(AUCrand_test))
	predictor_perf[drug,"meanAUCtrain"]=mean(AUCrand_train)
	predictor_perf[drug,"stdAUCtrain"]=std(t(AUCrand_train))
	predictor_perf[drug,"meanGI50"]=mean_cutoff
	
} #end drug loop


#Write summary table to file
write.table(predictor_perf, file=outfile, sep="\t", row.names=FALSE)
