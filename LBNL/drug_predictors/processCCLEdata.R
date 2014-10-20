#This script will take CCLE data files and create separate input files needed to run the drug predictor Rtoolbox
#To be used for downstream validation of drug predictors

#Load sample data
sample_info_file="/Users/ogriffit/Dropbox/drug_predictors/CCLE/CCLE_sample_info_file_2012-04-06.txt"
sample_info_data=read.table(file=sample_info_file, header=TRUE, sep="\t", as.is=c(1:13),row.names=1)

#Identify breast subset
breast_sample_info_data=sample_info_data[sample_info_data[,"Site.Primary"]=="breast",]
breast_samples=rownames(breast_sample_info_data)
breast_cell_line_names=breast_sample_info_data[,"Cell.line.primary.name"]

#Load expression data (Affymetrix U133+2 - pre-normalized with RMA, and summarized to gene symbol)
expression_data_file="/Users/ogriffit/Dropbox/drug_predictors/CCLE/CCLE_Expression_Entrez_2012-04-06.gct"
raw_expression_data=read.table(file=expression_data_file, header=TRUE, sep="\t", skip=2, as.is=c(1), row.names=1, comment.char="")
expression_gene_names=as.vector(raw_expression_data[,"Description"])

#Reduce expression data to single (most variable?) probeset for each gene
#It appears this data has already been summarized to gene level, however a few genes have probe id but not a symbol/description
#And, one gene symbol (TTL) appears for two probes
expression_gene_names[which(duplicated(expression_gene_names))]
apply(raw_expression_data[which(raw_expression_data[,"Description"]=="TTL"),2:length(colnames(raw_expression_data))],1,var)

#Exclude rows with symbol="NA" and choose just one probe for TTL
exclude=which(is.na(raw_expression_data[,"Description"]) | rownames(raw_expression_data)=="646982_at")
expression_data_filt=raw_expression_data[-exclude,]
expression_gene_names_filt=expression_gene_names[-exclude]

#Load copy number data (Affymetrix SNP6.0 - CBS segmentation, summarized at gene level)
CNV_data_file="/Users/ogriffit/Dropbox/drug_predictors/CCLE/CCLE_copynumber_byGene_2012-04-06.txt"
raw_CNV_data=read.table(file=CNV_data_file, header=TRUE, sep="\t", as.is=c(1:2), row.names=1, comment.char="")
CNV_gene_names=rownames(raw_CNV_data)

#Extract breast subset of data
#First determine which breast cell lines are present in both expression and copy number data
breast_samples2=breast_samples[which(breast_samples %in% colnames(raw_expression_data) & breast_samples %in% colnames(raw_CNV_data))]
breast_expression_data=expression_data_filt[,breast_samples2]
breast_CNV_data=raw_CNV_data[,breast_samples2]

#Set directory for breast specific files to be used in validation analysis
setwd("/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/TumorData")
#For each sample, create a separate expression and CNV file
#Also, create a mapping file with patient ids linked to expression and CNV file names
mapping=data.frame(CCLEsample=breast_samples2, Expfilename="NA", SNPfilename="NA", Methfilename="NA", row.names=breast_samples2, stringsAsFactors=FALSE)
for (i in 1:length(breast_samples2)){
	expdata=cbind(expression_gene_names_filt,breast_expression_data[,breast_samples2[i]])
	colnames(expdata)=c("Gene symbol",breast_samples2[i])
	expdata_file=paste("Affymetrix_CCLE_",breast_samples2[i],".txt",sep="")
	cnvdata=cbind(CNV_gene_names,breast_CNV_data[,breast_samples2[i]])
	colnames(cnvdata)=c("Gene symbol",breast_samples2[i])
	cnvdata_file=paste("SNP6_CCLE_",breast_samples2[i],".txt",sep="")
	write.table(expdata, file=expdata_file, quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)
	write.table(cnvdata, file=cnvdata_file, quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)
	mapping[breast_samples2[i],"Expfilename"]=expdata_file
	mapping[breast_samples2[i],"SNPfilename"]=cnvdata_file
}
setwd("/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/")
write.table(mapping, file="CCLE_affymetrix_SNP6_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
