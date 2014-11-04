datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/Affy"
#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/TMM"
setwd(datadir)

#Postscript/pdf files
#ps_file="LogFreq_dist_rand_vs_norm.ps"
#pdf(file="CommonModulesAffy_vs_random.pdf", width=7, height=7)
#pdf(file="CommonModules_AB2_vs_AB3plus_vs_DN_vs_random_AFFY.3.pdf", width=7, height=7)
#pdf(file="CommonModules_AB2_vs_AB3plus_vs_DN_vs_random_TMM.3.pdf", width=7, height=7)
pdf(file="CommonModules_AB2_vs_AB3plus_vs_DN_vs_random_AFFY_ClusterMeans.pdf", width=7, height=7)

#Get Coexpression data for genes defined by modules or random
#Affy (by gene pairs)
##data_rand = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_allgenes.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB2 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_doubletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
##data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPlusPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_DN = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Affy (by gene cluster)
data_rand = read.table("CommonModuleGeneClusters_allgenes.Affy.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_AB2 = read.table("AnnotationBased/clustermeans/CommonModuleGeneClusters_AnnotationBased_doubletPatterns.Affy.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_AB3 = read.table("AnnotationBased/clustermeans/CommonModuleGeneClusters_AnnotationBased_tripletPlusPatterns.Affy.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_DN = read.table("profileThresholded_250K/clustermeans/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#TMM (by gene pairs)
#data_rand = read.table("CommonModuleGeneClusters_allgenes.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB2 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_doubletPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPlusPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_DN = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Specify a maximum cluster size
maxclustersize=20
#maxclustersize=10
data_norm_AB3_undermax = data_norm_AB3[data_norm_AB3[,1]<=maxclustersize,][,2]
data_norm_AB2_undermax = data_norm_AB2[data_norm_AB2[,1]<=maxclustersize,][,2]
data_norm_DN_undermax = data_norm_DN[data_norm_DN[,1]<=maxclustersize,][,2]
data_rand_undermax = data_rand[data_rand[,1]<=maxclustersize,][,2]

#The following function can be used to get the scores and cumulative frequencies for the plot
cumfreq = function(x)
{
  x=sort(x)
  scorecounts=table(x)
  scores=names(scorecounts)
  cum_freqs=cumsum(scorecounts/sum(scorecounts))
  num_diff_scores=length(scores)
  result = matrix(0, nrow=num_diff_scores, ncol=2, dimnames=list(1:num_diff_scores, c("scores","cum_freqs")))
  result = cbind(scores,cum_freqs)
  return(result)
}

#Use cumfreq function to get scores and frequencies for each dataset to be plotted
results_AB2=cumfreq(data_norm_AB2_undermax)
results_AB3=cumfreq(data_norm_AB3_undermax)
results_DN=cumfreq(data_norm_DN_undermax)
results_rand=cumfreq(data_rand_undermax)


#Prepare to plot cum fractions for each on one plot
#plot(c(results_AB2[,"scores"],results_AB3[,"scores"],results_DN[,"scores"],results_rand[,"scores"]), c(results_AB2[,"cum_freqs"],results_AB3[,"cum_freqs"],results_DN[,"cum_freqs"],results_rand[,"cum_freqs"]), type="n", xlim=c(0,16), ylim=c(0.90,1), xlab="Coexpression score (TMM)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
plot(c(results_AB2[,"scores"],results_AB3[,"scores"],results_DN[,"scores"],results_rand[,"scores"]), c(results_AB2[,"cum_freqs"],results_AB3[,"cum_freqs"],results_DN[,"cum_freqs"],results_rand[,"cum_freqs"]), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Coexpression score (Pearson)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)

#Add points
points(results_AB2[,"scores"], results_AB2[,"cum_freqs"], pch=19, cex=0.25, col="orange")
points(results_AB3[,"scores"], results_AB3[,"cum_freqs"], pch=22, cex=0.25, col="darkgreen")
points(results_DN[,"scores"], results_DN[,"cum_freqs"], pch=23, cex=0.25, col="blue")
points(results_rand[,"scores"], results_rand[,"cum_freqs"], pch=25, cex=0.25, col="red")

#Plot lines if necessary
#lines(results_AB2[,"scores"], results_AB2[,"cum_freqs"], cex=0.25, col="orange")
#lines(results_AB3[,"scores"], results_AB3[,"cum_freqs"], cex=0.25, col="darkgreen")
#lines(results_DN[,"scores"], results_DN[,"cum_freqs"], cex=0.25, col="blue")
#lines(results_rand[,"scores"], results_rand[,"cum_freqs"], cex=0.25, col="red")


#Plot a legend
legend=c("Annotation Based (doublet)", "Annotation Based (triplet+)", "De Novo (triplet)", "Random")
legend(0.3,0.7, legend, pch=c(20,15,18,17), col=c("orange","darkgreen","blue","red"), cex=1.25)
#legend(5,0.94, legend, pch=c(20,15,18,17), col=c("orange","darkgreen","blue","red"), cex=1.25)

dev.off()

#Test for significance between distributions
#Problem is this tie warning an issue?
#ks.test(x=coexp_norm_AB3_sorted, y=coexp_rand_sorted, alternative = c("two.sided"), exact = NULL)
#ks.test(x=coexp_norm_DN_sorted, y=coexp_rand_sorted, alternative = c("two.sided"), exact = NULL)
