#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/SAGE"
datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/TMM"
setwd(datadir)

rownames=c(2-20)

### SAGE ###
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.SAGE_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.SAGE_by_clustersize.random100_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)

### cDNA ###
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.cDNA_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.cDNA_by_clustersize.random100_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)

### AFFY ###
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.Affy_by_clustersize.random100_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)

### TMM ###
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.random100_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=1, fill=T)

#With all individual coexpression values not averaged by clusters
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.3.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=rownames, fill=T)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.4.test.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=rownames, fill=T)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.4.test2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=rownames, fill=T)
data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.5.txt", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=rownames, fill=T)
#data_norm = read.table("temp", header=F, quote="", sep="\t", comment.char="", as.is=1, row.names=rownames, fill=T)

data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is, row.names=rownames, fill=T)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.random10000_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=rownames, row.names=1, fill=T)



num_cluster_sizes = length(data_norm[,1]) #Determines number of cluster sizes to loop through
cluster_sizes = dimnames(data_norm)[[1]] #Gets actual sizes for labels

#Create a vector to store stats in
cluster_stats = array(0, dimnames = list(cluster_sizes, c("mean", "std_dev", "rand_mean", "rand_std_dev")), dim=c(length(cluster_sizes),4))

for (i in 1:num_cluster_sizes){
cluster_stats[i,"mean"] = mean(as.numeric(data_norm[i,]),na.rm=T)
cluster_stats[i,"std_dev"] = sd(as.numeric(data_norm[i,]),na.rm=T)
cluster_stats[i,"rand_mean"] = mean(as.numeric(data_rand[i,]),na.rm=T)
cluster_stats[i,"rand_std_dev"] = sd(as.numeric(data_rand[i,]),na.rm=T)
}

#boxplot(t(data_norm))

plot(x=cluster_sizes,y=cluster_stats[,1], type="n",ylim=c(0,0.1),main='Mean TMM score for gene clusters sharing a common motif', xlab='Gene Cluster Size', ylab='Mean TMM score')
points(cluster_sizes,cluster_stats[,1],pch=19)
lines(cluster_sizes,cluster_stats[,1],lty=1)
points(cluster_sizes,cluster_stats[,3],pch=20)
lines(cluster_sizes,cluster_stats[,3],lty=2)

legend(.8,9.25,c('Gene clusters sharing a common motif','Random gene clusters'),lty=c(1,2))

headers=as.matrix(c("cluster_size","mean","std_dev","rand_mean","rand_std_dev"))
write (file="CommonMotif_vs_Coexpression_summary.txt", t(headers), ncolumns=5)
output = cbind(cluster_sizes,cluster_stats)
write (file="CommonMotif_vs_Coexpression_summary.txt", t(output), ncolumns=5, append=TRUE)
