setwd("/home/obig/Projects/sub_space_clustering/expO/TermOverRepresentation_all")

#datafile="GSE2109_gcrma_mapped_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster5_mindim15_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_gcrma_mapped_mincluster2_mindim50_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_random1_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_random1_TermOverRep_data_For_FisherExact.txt"

#Updated data for all terms
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_data_For_FisherExact.txt"
#datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_random1_TermOverRep_data_For_FisherExact.txt"
datafile="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_random1_TermOverRep_data_For_FisherExact.txt"

#Read in datafile
data=read.table(datafile, header=TRUE, sep="\t")

#Define function to calculate fisher exact statistic for each cluster/term combination
cont_table_fun=function(x){
cont_table=matrix(x,nr = 2,dimnames=list(c("interm", "notinterm"),c("inclust", "notclust")))
fisher_test_result=fisher.test(cont_table, alternative = "greater")
result=matrix(c(fisher_test_result$p.value, fisher_test_result$estimate[[1]]), nr=1)
return(result)
}

#Apply the fisher function to all rows of the data for the appropriate columns
fisher_results=t(apply(data[,c(3,4,5,6)], 1, cont_table_fun))

#Give some column names to the output
colnames(fisher_results)=c("pvalue", "odds_ratio")

#reduce number of decimal places
#fisher_results_formatted=format(fisher_results, digits=5, scientific=FALSE)

#Combine fisher exact results with initial data
final_results=cbind(data,fisher_results)

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
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_TermOverRep_FisherExact_results.txt")
#write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster5_mindim15_random1_TermOverRep_FisherExact_results.txt")
write.table(final_results_ordered, row.names=FALSE, sep="\t", quote=FALSE, file="GSE2109_1026exps_16Aug2006_gcrma_ENSG_mappable2uniprot_AND_ENSG_overlap_mincluster2_mindim50_random1_TermOverRep_FisherExact_results.txt")
