library(randomForest)

#drug="BEZ235"
#drug="17-AAG"
#drug="BIBW2992"
#drug="GSK_Tykerb"
drug="GSK1838705A (IGF1R)"

RFmodelfile=paste("C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/independent/rnaseq/RFmodels/RF_model_RS_",drug,".Rdata",sep="")
RFoptmodelfile=paste("C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/independent/rnaseq/RFmodels/RF_model_opt_RS_",drug,".Rdata",sep="")
outfile=paste(drug,"_variables.txt",sep="")
outfile2=paste(drug,"_variables_diff.txt",sep="")


load(file=RFmodelfile)
load(file=RFoptmodelfile)
model_variables=names(rf_model$importance[,1])
variable_importances=rf_model$importance[,4]
opt_model_variables=names(rf_model_opt$importance[,1])
opt_variable_importances=rf_model_opt$importance[,4]

setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/splicing_importance/")

model_variables_split=strsplit(x=model_variables, split="__")
getPrefix=function(x){y=x[[1]][1]; return(y)}
prefixes=sapply(model_variables_split,getPrefix, simplify = TRUE)

#Make predictor names without prefix, drop leading ID for RNAseq Names (only type to have 3 components)
getName=function(x){y=x[2:length(x)]; return(y)} 

Names=sapply(model_variables_split,getName,simplify=FALSE)
pasteNames=function(x){
 if (length(x)==3){
  a=x[1]
  b=x[2]
  c=x[3]
 }else if (length(x)==2){
  a="NA"
  b=x[1]
  c=x[2]
 }else{
  a="NA"
  b="NA"
  c=x[1]
 }
y=paste(a,b,c,sep=",")
return(c(a,b,c))
}
CleanNames_list=sapply(Names,pasteNames,simplify=FALSE)  
CleanNames_matrix=do.call(rbind, CleanNames_list)

model_variables_out=cbind(prefixes,model_variables,CleanNames_matrix,variable_importances)
colnames(model_variables_out)=c("prefixes","model_variables","feature","feature2","gene","variable_importances")

write.table(model_variables_out, file=outfile, sep="\t", col.names=TRUE, row.names=FALSE)

#Get list of unique genes
genes=unique(model_variables_out[,"gene"])

#Create dataframe to store results
results=matrix(data=NA, nrow=length(genes), ncol=5, dimnames=list(genes, c("gene","gene_VI","best_feature","best_feature_VI","diff_VI")))

for (i in 1:length(genes)){
 #Extract data for gene
 gene=genes[i]
 print(gene)
 gene_data=model_variables_out[which(model_variables_out[,"gene"]==gene),,drop=FALSE]

 #Find best feature and compare to gene
 if(gene %in% gene_data[,"feature2"]){
  gene_var_imp=as.numeric(gene_data[which(gene_data[,"feature2"] %in% gene),"variable_importances"])
 }else{
  gene_var_imp="NA"
 }

 best_var=gene_data[order(as.numeric(gene_data[,"variable_importances"]), decreasing=TRUE)[1],"model_variables"]
 best_var_imp=as.numeric(gene_data[order(as.numeric(gene_data[,"variable_importances"]), decreasing=TRUE)[1],"variable_importances"])
 if(gene_var_imp=="NA"){
  var_imp_diff=best_var_imp
 }else{
  var_imp_diff=best_var_imp-gene_var_imp
 }

 #Add to results
 results[gene,"gene"]=gene
 results[gene,"gene_VI"]=gene_var_imp
 results[gene,"best_feature"]=best_var
 results[gene,"best_feature_VI"]=best_var_imp
 results[gene,"diff_VI"]=var_imp_diff
}

write.table(results, file=outfile2, sep="\t", col.names=TRUE, row.names=FALSE)




