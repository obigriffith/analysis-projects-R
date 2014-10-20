setwd("C:/Documents and Settings/obig/My Documents/Projects/subspace_clustering/KiWi_vs_cisRED")

kiwi_data_file="kiwi.dat"
random_data_file="random_new.dat"
#random_data_file="random.dat"

#Format for kiwi data is as follows: clusterId size dim aveScore
kiwi_data=read.table(file=kiwi_data_file, header=FALSE, sep="\t", row.names=1)

#Format is as follows: random_list_name size aveScore
#Note that only size was variable, the data is for all dimensions
random_data=read.table(file=random_data_file, header=FALSE, sep="\t", row.names=1)

#Get dimension numbers for kiwi data
kiwi_dims=sort(unique(kiwi_data[,2]))

#Get cluster sizes for random data
random_cluster_sizes=sort(unique(random_data[,1]))

#Get all combinations of dims and cluster sizes in data
kiwi_data_ordered=kiwi_data[order(kiwi_data[,2], kiwi_data[,1]),]
combos=unique(kiwi_data_ordered[,1:2])
combo_names=paste(combos[,1], combos[,2],sep="_")

#Create array to store results
kiwi_results = array(0, dimnames = list(combo_names, c("dims","cluster_size","num_clusters","mean_score","sd_score")), dim=c(length(combo_names),5))

#Summarize mean and sd of score for kiwi data for each dimension and clustersize
#Go through each dimension size
n=0
for (i in kiwi_dims){
 #Get KiWi with dimension=i
 kiwi_data_subset=kiwi_data[kiwi_data[,2]==i,]
 kiwi_cluster_sizes=sort(unique(kiwi_data_subset[,1]))
 #Go through each cluster size
 for (j in kiwi_cluster_sizes){
  n=n+1
  #Get scores for just this dimension and cluster size
  kiwi_score_subset=kiwi_data_subset[kiwi_data_subset[,1]==j,3]
  kiwi_sd_score=sd(kiwi_score_subset)
  kiwi_mean_score=mean(kiwi_score_subset)
  kiwi_results[n,"dims"]=i
  kiwi_results[n,"cluster_size"]=j
  kiwi_results[n,"num_clusters"]=length(kiwi_score_subset)
  kiwi_results[n,"mean_score"]=kiwi_mean_score
  kiwi_results[n,"sd_score"]=kiwi_sd_score
 }
}
write.table (kiwi_results, sep="\t", file="kiwi_results.txt")


#Summarize mean and sd of score for kiwi data for each clustersize (all dimensions)
kiwi_cluster_sizes=sort(unique(kiwi_data[,1]))
#Create array to store results for kiwi data
kiwi_alldim_results = array(0, dimnames = list(kiwi_cluster_sizes, c("num_clusters","mean_score", "sd_score")), dim=c(length(kiwi_cluster_sizes),3))
n=0
for (k in kiwi_cluster_sizes){
n=n+1
 #Get scores for just this dimension and cluster size
  kiwi_score_subset=kiwi_data[kiwi_data[,1]==k,3]
  kiwi_sd_score=sd(kiwi_score_subset)
  kiwi_mean_score=mean(kiwi_score_subset)
  kiwi_alldim_results[n,"num_clusters"]=length(kiwi_score_subset)
  kiwi_alldim_results[n,"mean_score"]=kiwi_mean_score
  kiwi_alldim_results[n,"sd_score"]=kiwi_sd_score
}
write.table (kiwi_alldim_results, sep="\t", file="kiwi_alldim_results.txt")

#Summarize mean and sd of score for random data for each clustersize (all dimensions)
#Create array to store results for kiwi data
random_alldim_results = array(0, dimnames = list(random_cluster_sizes, c("num_clusters","mean_score", "sd_score")), dim=c(length(random_cluster_sizes),3))
n=0
for (l in random_cluster_sizes){
n=n+1
 random_score_subset=random_data[random_data[,1]==l,2]
 random_sd_score=sd(random_score_subset)
 random_mean_score=mean(random_score_subset)
 random_alldim_results[n,"num_clusters"]=length(random_score_subset)
 random_alldim_results[n,"mean_score"]=random_mean_score
 random_alldim_results[n,"sd_score"]=random_sd_score
}
write.table (random_alldim_results, sep="\t", file="random_alldim_results.txt")


#Overall stats
kiwi_all_scores=kiwi_data[,3]
random_all_scores=random_data[,2]
mean(random_all_scores)
mean(kiwi_all_scores)
wilcox.test(x=kiwi_all_scores, y=random_all_scores, alternative ="two.sided", paired = FALSE, exact = NULL, correct = TRUE)

