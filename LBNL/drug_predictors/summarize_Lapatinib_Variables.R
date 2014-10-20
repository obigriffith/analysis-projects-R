library(randomForest)
load(file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/core_alldata_imputed/RFmodels/RF_model_opt_GSK_Tykerb")
model_variables=names(rf_model_opt$importance[,1])
variable_importances=rf_model_opt$importance[,4]

setwd("C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/allDrugs/core_alldata_imputed/")

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

#Define set of HER2Amp Genes
HER2Amp_Genes=c("FBXL20","PPARBP","MED1","CRKRS","NEUROD2","PPP1R1B","STARD3","TCAP","PNMT","PERLD1","PGAP3","ERBB2","HER2p1248","C17orf37","GRB7","IKZF3","ZNFN1A3","ZPBP2","GSDMB","ORMDL3","GSDMA","PSMD3","CSF3","THRAP4","DRIP100","MED24","THRA","NR1D1","MSL1","CASC3","RAPGEFL1")
HER2Amp_Genes_core=c("NEUROD2","PPP1R1B","STARD3","TCAP","PNMT","PERLD1","PGAP3","ERBB2","HER2p1248","C17orf37","GRB7","IKZF3","ZNFN1A3")


model_variables_out=cbind(model_variables_out,model_variables_out[,5]%in%HER2Amp_Genes_core,model_variables_out[,5]%in%HER2Amp_Genes)
colnames(model_variables_out)=c("data_type","model_variable","RS_feature","feature","Gene","VarImp","HER2ampCore","HER2ampExpanded")

write.table(model_variables_out, file="RF_model_opt_GSK_Tykerb_variables.txt", sep="\t", col.names=TRUE, row.names=FALSE)




#Get ERBB2 Amp genes
### exp set - flanking (500,000bp up-stream) ###
#"FBXL20"
#"PPARBP" #aka MED1
#"CRKRS"
### exp set ### - core amplicon
#"NEUROD2"
#"PPP1R1B"
#"STARD3"
#"TCAP"
#"PNMT"
#"PERLD1" #AKA PGAP3
#"ERBB2"
#"C17orf37"
#"GRB7"
#"IKZF3" #AKA ZNFN1A3
### exp set ### - core amplicon
#"ZPBP2"
#"GSDMB"
#"ORMDL3"
#"GSDMA"
#"PSMD3"
#"CSF3"
#"THRAP4" #aka DRIP100, MED24
#"THRA"
#"NR1D1"
#"MSL1"
#"CASC3"
#"RAPGEFL1"
### exp set ### - flanking (500,000bp down-stream) ###

