#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/SAGE"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/cDNA"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/Affy"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/TMM"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/TMM/cisred_1_2e"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/Affy/cisred_1_2e"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/cDNA/cisred_1_2e"
#datadir = "~/Projects/cisRED_coexpression/commonMotifAnalysis/SAGE/cisred_1_2e"

#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/SAGE"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cDNA"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/Affy"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/TMM"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/TMM/cisred_1_2e"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/TMM"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/Affy/profileThresholded_250K/"
datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/Affy/AnnotationBased"

setwd(datadir)

#Results file
#outfile="CommonMotifGeneClusters_results_SAGE_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_results_cDNA_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_results_Affy_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_results_TMM_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_1_2e_results_TMM_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_1_2e_results_Affy_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_1_2e_results_cDNA_vs_random10000pairs_20.txt"
#outfile="CommonMotifGeneClusters_1_2e_results_SAGE_vs_random10000pairs_20.txt"

#outfile="CommonModuleGeneClusters_results_TMM_vs_random10000pairs_20.2.txt"
#outfile="CommonModuleGeneClusters_results_Affy_vs_random10000pairs_20.2.txt"
#outfile="CommonModuleGeneClusters_results_cDNA_vs_random10000pairs_20.2.txt"
#outfile="CommonModuleGeneClusters_results_SAGE_vs_random10000pairs_20.2.txt"
#outfile="CommonModuleGeneClusters_1_2e_results_TMM_vs_random10000pairs_20.2.txt"
#outfile="CommonModuleGeneClusters_2_results_TMM_vs_random10000pairs_20.txt"
outfile="CommonModuleGeneClusters_2_AnnotationBased_results_Affy_vs_random10000pairs_25.txt"

#Postscript file
ps_file="LogFreq_dist_rand_vs_norm.ps"

#Common motif analysis
#SAGE
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.SAGE_by_clustersize.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.SAGE_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2e.SAGE_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2e.SAGE_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#cDNA
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.cDNA_by_clustersize.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.cDNA_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2e.cDNA_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2e.cDNA_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#Affy
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.Affy_by_clustersize.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2e.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2e.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#TMM
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.6.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.random10000pairs_20.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonMotifGeneClusters_cisred_1_2e.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonMotifGeneClusters_cisred_1_2e.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Common Module Analysis #can use the same random data as these are just randomly create gene clusters from cisRED genes
#TMM
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.6.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.TMM_by_clustersize.txt.2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("~/Projects/cisRED_coexpression/commonMotifAnalysis/TMM/CommonMotifGeneClusters_cisred_1_2a.TMM_by_clustersize.random10000pairs_20.2.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2e.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_cisred_1_2e.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Affy
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.6.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("~/Projects/cisRED_coexpression/commonMotifAnalysis/Affy/CommonMotifGeneClusters_cisred_1_2a.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.random10000pairs_30.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm = read.table("CommonModuleGeneClusters_AnnotationBased_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_rand = read.table("CommonModuleGeneClusters_AnnotationBased_tripletPatterns.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#cDNA
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.cDNA_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.6.cDNA_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("~/Projects/cisRED_coexpression/commonMotifAnalysis/cDNA/CommonMotifGeneClusters_cisred_1_2a.cDNA_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#SAGE
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.SAGE_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm = read.table("CommonModuleGeneClusters_cisred_1_2a.6.SAGE_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("~/Projects/cisRED_coexpression/commonMotifAnalysis/SAGE/CommonMotifGeneClusters_cisred_1_2a.SAGE_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Set cluster sizes
clustersizes=c(5,6,7,8,10,11,12,15,24)
#clustersizes=c(2:25)
#clustersizes=c(2:20)
#clustersizes=c(2:6)

#Create a vector to store stats in
cluster_stats = array(0, dimnames = list(clustersizes, c("n", "mean", "std_dev", "std_err", "n_rand", "rand_mean", "rand_std_dev", "rand_std_err", "ttest")), dim=c(length(clustersizes),9))

num_cluster_sizes = length(clustersizes)

for (i in 1:num_cluster_sizes){
clustersize=clustersizes[i]
data_norm_for_clustersize = data_norm[data_norm[,1]==clustersize,]
data_rand_for_clustersize = data_rand[data_rand[,1]==clustersize,]
cluster_stats[i,"n"] = length(data_norm_for_clustersize[,2])
cluster_stats[i,"mean"] = mean(data_norm_for_clustersize[,2],na.rm=T)
cluster_stats[i,"std_dev"] = sd(data_norm_for_clustersize[,2],na.rm=T)
cluster_stats[i,"std_err"] = cluster_stats[i,"std_dev"]/sqrt(cluster_stats[i,"n"] - 1) #SE = sd/SQRT(n - 1)
cluster_stats[i,"n_rand"] = length(data_rand_for_clustersize[,2])
cluster_stats[i,"rand_mean"] = mean(data_rand_for_clustersize[,2],na.rm=T)
cluster_stats[i,"rand_std_dev"] = sd(data_rand_for_clustersize[,2],na.rm=T)
cluster_stats[i,"rand_std_err"] = cluster_stats[i,"rand_std_dev"]/sqrt(cluster_stats[i,"n_rand"] - 1) #SE = sd/SQRT(n - 1)
norm_vs_rand_ttest=t.test(x=data_norm_for_clustersize[,2],y=data_rand_for_clustersize[,2], alternative=c("two.sided"), mu=0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
cluster_stats[i,"ttest"] = norm_vs_rand_ttest$p.value
}


headers=as.matrix(c("n", "mean", "std_dev", "std_err", "n_rand", "rand_mean", "rand_std_dev", "rand_std_err", "ttest"))
write (file=outfile, t(headers), ncolumns=10)
output = cbind(clustersizes,cluster_stats)
write (file=outfile, t(output), ncolumns=10, append=T)

#To plot distribution of normal vs random data for all clustersizes
postscript(ps_file, pointsize=1, width=9.5, height=7)
hist_df = hist(data_norm[,2], nclass=1000, plot=FALSE)
hist_df_rand = hist(data_rand[,2], nclass=1000, plot=FALSE)
par(mai=c(0.5,0.5,0.5,0.5))
plot(c(hist_df$mids,hist_df_rand$mids), c(hist_df$density,hist_df_rand$density), type="n", log="y", xlim=c(0,1), xlab="Coexpression score", ylab="log frequency", main="frequency vs coexpression score for gene pairs with a common module",cex.lab=2.5, cex.axis=2.5, cex.main=2)
#plot(c(hist_df$mids,hist_df_rand$mids), c(hist_df$density,hist_df_rand$density), type="n", log="y", xlim=c(0,10), xlab="Coexpression score", ylab="log frequency", main="frequency vs coexpression score for gene pairs with a common module",cex.lab=2.5, cex.axis=2.5, cex.main=2)
points(hist_df$mids, hist_df$density, pch=19, cex=2.5)
points(hist_df_rand$mids, hist_df_rand$density, pch=23, cex=2.5)
legend=c("norm","rand")
#legend(8,8, legend, pch=c(19,23),cex=2.5)
legend(0.8,100, legend, pch=c(19,23),cex=2.5)
dev.off()


#To plot distribution of normal vs random data for any clustersize
#clustersize=2
#data_norm_for_clustersize = data_norm[data_norm[,1]==clustersize,]
#data_rand_for_clustersize = data_rand[data_rand[,1]==clustersize,]
#hist_df = hist(data_norm_for_clustersize[,2], nclass=1000, plot=FALSE)
#hist_df_rand = hist(data_rand_for_clustersize[,2], nclass=1000, plot=FALSE)
#plot(c(hist_df$mids,hist_df_rand$mids), c(hist_df$density,hist_df_rand$density), type="n", log="y", xlim=c(0,6), xlab="Coexpression score", ylab="log frequency", main="frequency vs coexpression score (TMM) for gene pairs sharing a common motif", cex=2)
#points(hist_df$mids, hist_df$density, pch=19)
#points(hist_df_rand$mids, hist_df_rand$density, pch=23)
#legend=c("norm","rand")
#legend(0,1, legend, pch=c(19,23))
