library(randomForest)
library(mclust)

outdir="C:/Users/Obi/Documents/Projects/BCCL/drug_response_predictors/"
outfile="predictor_performance_genes.txt"
#outfile="predictor_performance_junctions.txt"
setwd(outdir)

#Parameters
pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 2 #Minimum coefficient of variation (0.7 recommended?)
cov_max = 10 #Maximum cov

#Libraries that have drug response, RNAseq and mutation data
complete_libs=c("BT20","BT474","BT549","CAMA1","HCC1143","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC202","HCC3153","HCC38","HCC70","HS578T","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SUM1315MO2","SUM149PT","SUM159PT","T47D","UACC812","UACC893","ZR751","ZR7530")

#Libraries that have drug response and RNAseq data
RNAseq_libs=c("X184A1","X184B5","X600MPE","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC3153","HCC38","HCC70","HS578T","LY2","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SKBR3","SUM1315","SUM149PT","SUM159PT","SUM225CWN","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")


#RNAseq data (Choose one feature type)
#Gene
datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_GeneExpression_v53.txt"

#Transcript
#datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_TranscriptExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_TranscriptExpression_v53.txt"

#Junction
#datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Matrix_KnownJunctionExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/57libs/matrix/Expressed_KnownJunctionExpression_v53.txt"


raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:62]
raw_exp_status=raw_exp_status_import[,5:62]

#Fix misnamed library
colnames(raw_data)[which(colnames(raw_data)=="MDAMB13v1")]="MDAMB134VI"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="MDAMB13v1")]="MDAMB134VI"

#Retrieve data for only libraries with both RNAseq and drug response data
data=raw_data[,RNAseq_libs]
feat_data=raw_feat_data
exp_status=raw_exp_status[,RNAseq_libs]
libs=colnames(data)

#Define a percent expressed function and filter out features with less than minimum
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}
pe_data=apply(exp_status, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data_filt=data[passed_pe,]
feat_data_filt=feat_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
  return(cov)
}

cov_data=apply(data_filt, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data_filt=data_filt[passed_cov,]
feat_data_filt=feat_data_filt[passed_cov,]
rownames(data_filt)=feat_data_filt[,"FID"]

#drug response data
drugdatafile="C:/Users/Obi/Documents/Projects/BCCL/drug_data/LBNL_Gray_BCCL_GI50_10Feb_final_encoded.txt"
raw_drugdata_import=read.table(drugdatafile, header = TRUE, na.strings = "NA", sep="\t")

drug_data=raw_drugdata_import[,2:75]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Retrieve data for only libraries with both RNAseq and drug response data
drug_data_filt=drug_data[RNAseq_libs,]

#Transform data to -log values
drug_data_filt_trans=-log(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#Apply COV filter???
#cov_data2=apply(drug_data_filt_trans, 2, cov_fun)


#For each drug, find expression features associated with response and build predictor for drug response
#First, create dataframe to hold results
predictor_perf = data.frame(cbind(drugs,acc=NA, err=NA), stringsAsFactors=FALSE)

for (i in 1:length(drugs)){
  #Divide cell lines into responders and non-responders for drug
  drug=drugs[i]
  x=as.numeric(drug_data_filt_trans[,i])
  x_NA=which(is.na(x))
  x_ordered=order(x,na.last=NA)
  cutoff=round(length(x_ordered)/2)
  resistants=x_ordered[1:cutoff]
  sensitives=x_ordered[(cutoff+1):length(x_ordered)]
  response_class=vector(length=length(x))
  response_class[x_NA]=NA
  response_class[sensitives]="sensitive"
  response_class[resistants]="resistant"

  #Exclude libs where response_class=NA
  nonNA=which(!is.na(response_class))
  data_filt_nonNA=data_filt[,nonNA]
  response_class_nonNA=response_class[nonNA]

  #Make sure there are at least 10 libs classified as sensitive and resistant
  if (length(which(response_class_nonNA=="sensitive")) >= 10 & length(which(response_class_nonNA=="resistant")) >= 10){

    #Create RF predictor for response_class using RNAseq data(data_filt)
    rf_output=randomForest(x=t(data_filt_nonNA), y=as.factor(response_class_nonNA), importance = TRUE, ntree = 1001, proximity=TRUE)

    #Determine performance statistics
    overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
    overall_accuracy=100-overall_error
    predictor_perf[i,"err"]=overall_error
    predictor_perf[i,"acc"]=overall_accuracy
  }
}

write.table(predictor_perf, file=outfile, sep="\t", row.names=FALSE)






#Unused code
#Attempt to define cutoffs by mixed-model clustering of drug responses
#Force mclust to break into 3 clusters 
#x=as.numeric(drug_data_filt_trans[,i])
#x=x[which(!is.na(x))]
#mclust_response=Mclust(x, G=3, na.rm=TRUE)
#summary(mclust_response, x) #Gives you list of values returned
#classification=mclust_response$classification

#Choose cutoffs - use max and min of "middle" cluster (#2) to break into three components
#mm_cutoff1=min(x[classification==2])
#mm_cutoff2=max(x[classification==2])



