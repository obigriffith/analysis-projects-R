#Load model to extract features from
RFmodeldir="/Users/ogriffit/Archive/LBNL_linux/Projects/drug_predictors/RFpredictors/allDrugs/independent/gene_level/rnaseq/RFmodels/"
#models=list.files(path=RFmodeldir, include.dirs=FALSE, pattern="_opt") # Get opt model files only
models=list.files(path=RFmodeldir, include.dirs=FALSE) # Get all model files
model_count=length(models)

signaturedir="/Users/ogriffit/Archive/LBNL_linux/Projects/drug_predictors/RFpredictors/allDrugs/independent/gene_level/rnaseq/RFmodels/signatures"
setwd(signaturedir)

for(i in 1:model_count){
	model=models[i]
	print(paste("processing",model,sep=" "))
	RFmodelpath=paste(RFmodeldir,model,sep="")
	load(RFmodelpath)
	gene_names=rownames(rf_model_opt$importance)
	#gene_names=gsub("__",",",gene_names)
	sig_genes=cbind(gene_names,rf_model_opt$importance[,"MeanDecreaseGini"],rf_model_opt$importance[,"MeanDecreaseAccuracy"])
	colnames(sig_genes)=c("Feature","VarImp","VarImpSD")
	outfile=gsub("Rdata","txt", model)
	write.table(sig_genes, file=outfile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

