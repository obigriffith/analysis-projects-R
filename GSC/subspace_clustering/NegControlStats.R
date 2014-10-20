setwd("C:/Documents and Settings/obig/My Documents/Projects/subspace_clustering/StanfordPromoters/NegControlAnalysis")

kiwi_data_file="cooper.expander.shortNames.probeId.20060328_minclust2_mindim6_NegControlSummary.3.txt"
random_data_file="cooper.expander.shortNames.probeId.20060328_minclust2_mindim6_random.1.NegControlSummary.3.txt"

#Format for data is as follows: clusterId num_dims num_genes num_negs
kiwi_data=read.table(file=kiwi_data_file, header=FALSE, sep="\t", row.names=1)
random_data=read.table(file=random_data_file, header=FALSE, sep="\t", row.names=1)

#Get dimension numbers and cluster sizes for kiwi data (in this case, random should be exactly the same)
kiwi_dims=sort(unique(kiwi_data[,1]))
kiwi_sizes=sort(unique(kiwi_data[,2]))

#Get all combinations of dims and cluster sizes in data
kiwi_data_ordered=kiwi_data[order(kiwi_data[,1], kiwi_data[,2]),]
combos=unique(kiwi_data_ordered[,1:2])
combo_names=paste(combos[,1], combos[,2],sep="_")

#Create array to store results
kiwi_results = array(0, dimnames = list(combo_names, c("dims","cluster_size","num_clusters","mean_score","sd_score")),dim=c(length(combo_names),5))

#Summarize mean and sd of score for kiwi data for each dimension and clustersize
#Go through each dimension size
n=0
for (i in kiwi_dims){
 #Get KiWi with dimension=i
 kiwi_data_subset=kiwi_data[kiwi_data[,1]==i,]
 kiwi_cluster_sizes=sort(unique(kiwi_data_subset[,2]))
 #Go through each cluster size
 for (j in kiwi_cluster_sizes){
  n=n+1
  #Get scores for just this dimension and cluster size
  kiwi_score_subset=kiwi_data_subset[kiwi_data_subset[,2]==j,]
  #Get fraction/percent of genes that are negative controls per cluster
  kiwi_score_subset_fractions=kiwi_score_subset[,3]/kiwi_score_subset[,2]
  kiwi_sd_score=sd(kiwi_score_subset_fractions)
  kiwi_mean_score=mean(kiwi_score_subset_fractions)
  kiwi_results[n,"dims"]=i
  kiwi_results[n,"cluster_size"]=j
  kiwi_results[n,"num_clusters"]=length(kiwi_score_subset_fractions)
  kiwi_results[n,"mean_score"]=kiwi_mean_score
  kiwi_results[n,"sd_score"]=kiwi_sd_score
 }
}

write.table (kiwi_results, sep="\t", file="kiwi_results.txt")


#Summarize mean and sd of score for kiwi data for each clustersize (all dimensions)
kiwi_cluster_sizes=sort(unique(kiwi_data[,2]))
#Create array to store results for kiwi data
kiwi_alldim_results = array(0, dimnames = list(kiwi_cluster_sizes, c("num_clusters","mean_score", "sd_score")), dim=c(length(kiwi_cluster_sizes),3))
n=0
for (k in kiwi_cluster_sizes){
n=n+1
 #Get scores for just this dimension and cluster size
  kiwi_score_subset=kiwi_data[kiwi_data[,2]==k,]
  kiwi_score_subset_fractions=kiwi_score_subset[,3]/kiwi_score_subset[,2]
  kiwi_sd_score=sd(kiwi_score_subset_fractions)
  kiwi_mean_score=mean(kiwi_score_subset_fractions)
  kiwi_alldim_results[n,"num_clusters"]=length(kiwi_score_subset_fractions)
  kiwi_alldim_results[n,"mean_score"]=kiwi_mean_score
  kiwi_alldim_results[n,"sd_score"]=kiwi_sd_score
}
write.table (kiwi_alldim_results, sep="\t", file="kiwi_alldim_results.txt")

#Summarize mean and sd of score for random data for each clustersize (all dimensions)
#Create array to store results for kiwi data
random_cluster_sizes=sort(unique(random_data[,2]))
random_alldim_results = array(0, dimnames = list(random_cluster_sizes, c("num_clusters","mean_score", "sd_score")), dim=c(length(random_cluster_sizes),3))
n=0
for (l in random_cluster_sizes){
n=n+1
 random_score_subset=random_data[random_data[,2]==l,]
 random_score_subset_fractions=random_score_subset[,3]/random_score_subset[,2]
 random_sd_score=sd(random_score_subset_fractions)
 random_mean_score=mean(random_score_subset_fractions)
 random_alldim_results[n,"num_clusters"]=length(random_score_subset_fractions)
 random_alldim_results[n,"mean_score"]=random_mean_score
 random_alldim_results[n,"sd_score"]=random_sd_score
}
write.table (random_alldim_results, sep="\t", file="random_alldim_results.txt")


#Test whether distribution of negative fractions is different between KiWi and Random
wilcox.test(x=kiwi_data[,3]/kiwi_data[,2], y=random_data[,3]/random_data[,2], paired=FALSE, alternative="two.sided")
#W = 21013911365, p-value < 2.2e-16 (Mann-Whitney U-test)