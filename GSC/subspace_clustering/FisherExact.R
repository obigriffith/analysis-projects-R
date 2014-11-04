#setwd("/home/obig/Projects/sub_space_clustering/expO/TermOverRepresentation/")
#setwd("/home/obig/Projects/sub_space_clustering/expO/TermOverRepresentation_all/")
setwd("/home/obig/Projects/sub_space_clustering/expO/TermOverRepresentation_test/")

#datafile="GSE2109_gcrma_mapped_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster5_mindim15_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster2_mindim50_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_random1_TermOverRep_data_For_FisherExact.txt"

#Updated data for all terms
datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"

data=read.table(datafile, header=TRUE, sep="\t")

cluster_names=data[,1]
terms=data[,2]

#Create array to store results
pvalue_results = array(0, dimnames = list(cluster_names, c("pvalue", "odds_ratio")), dim=c(length(cluster_names),2))

for (i in 1:length(cluster_names)){
cont_table=matrix(c(data[i,3], data[i,4], data[i,5], data[i,6]),nr = 2,dimnames=list(c("interm", "notinterm"),c("inclust", "notclust")))
fisher_test_result=fisher.test(cont_table, alternative = "greater")
#Enter results into array for print out
pvalue_results[i,"pvalue"] =  fisher_test_result$p.value
pvalue_results[i,"odds_ratio"] = fisher_test_result$estimate
}

#reduce number of decimal places
#pvalue_results_formatted=format(pvalue_results, digits=5, scientific=FALSE)

#Combine fisher exact results with initial data
final_results=cbind(data,pvalue_results)

#Order by pvalue
final_results_ordered=final_results[order(final_results$pvalue),]

#Write results to file
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_gcrma_mapped_mincluster5_mindim15_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_gcrma_mapped_mincluster5_mindim15_random1_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_gcrma_mapped_mincluster2_mindim50_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_gcrma_mapped_mincluster2_mindim50_random1_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_random1_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_random1_TermOverRep_FisherExact_results.txt")

#Updated results for all terms
write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_FisherExact_results.txt")
