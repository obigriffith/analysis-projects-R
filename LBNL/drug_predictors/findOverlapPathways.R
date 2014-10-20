

datadir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ClueGO/"
setwd(datadir)

pathway_data=read.table(file="CompleteListPathways_vs_Drugs_noRPPA.txt", header=TRUE, as.is=c(1:3), sep="\t")

#Summarize counts for each term
termcounts=table(pathway_data[,2])
unique_terms=names(termcounts)

#Get corresponding decription for each term
term2Desc=unique(pathway_data[,2:3])
rownames(term2Desc)=term2Desc[,1]
unique_descs=term2Desc[unique_terms,2]

#Get list of drugs for each term
getdrugs=function(x){
 drugs=paste(pathway_data[which(pathway_data[,"GOID"]==x),"Drug"],collapse=",")
 return(drugs)
}
drug_lists=sapply(unique_terms, getdrugs, simplify=TRUE)
result=cbind(as.vector(unique_terms), as.vector(termcounts), as.vector(unique_descs), as.vector(drug_lists))
colnames(result)=c("Term","Count","Description","Drugs")

write.table(result, file="PathwaysDrugCounts_noRPPA.txt", sep="\t", quote=TRUE, row.names=FALSE)