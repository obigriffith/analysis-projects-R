setwd("C:/Documents and Settings/obig/My Documents/Projects/subspace_clustering/ClusterRepeatAnalysis")

kiwi_data_file="GSE2109_1026exps_16Aug2006_gcrma_probe_mappable2uniprot_AND_ENSG_overlap_repeatanalysis.2.txt"

#Format for data is as follows: cluster_number num_genes num_repeats
kiwi_data=read.table(file=kiwi_data_file, header=FALSE, sep="\t", row.names=1)

#Get dimension numbers and cluster sizes for kiwi data
kiwi_sizes=sort(unique(kiwi_data[,1]))

#Create array to store results
kiwi_results = array(0, dimnames = list(kiwi_sizes, c("num_clusters","mean_score","sd_score")),dim=c(length(kiwi_sizes),3))

#Summarize mean and sd of score for kiwi data for each clustersize (all dimensions)
n=0
for (k in kiwi_sizes){
n=n+1
 #Get scores for just this dimension and cluster size
  kiwi_score_subset=kiwi_data[kiwi_data[,1]==k,]
  kiwi_score_subset_scores=kiwi_score_subset[,2]
  kiwi_sd_score=sd(kiwi_score_subset_scores)
  kiwi_mean_score=mean(kiwi_score_subset_scores)
  kiwi_results[n,"num_clusters"]=length(kiwi_score_subset_scores)
  kiwi_results[n,"mean_score"]=kiwi_mean_score
  kiwi_results[n,"sd_score"]=kiwi_sd_score
}

write.table (kiwi_results, sep="\t", file="kiwi_results.txt")
