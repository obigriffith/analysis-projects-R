datadir = "~/Projects/cisRED_coexpression/antisense_RNA/"
setwd(datadir)

#Postscript/pdf files
#pdf(file="antisense_RNA_1_vs_2_vs_3_vs_rand_Affy.pdf", width=7, height=7)
#pdf(file="antisense_RNA_1_vs_2_vs_3_vs_rand_Affy.2.pdf", width=7, height=7)
#pdf(file="antisense_RNA_1_vs_2_vs_3_vs_rand_Affy.3.pdf", width=7, height=7)
pdf(file="antisense_RNA_1_vs_2_vs_3_vs_rand_Affy_nonzero.1.pdf", width=7, height=7)

#Get Coexpression data for genes pairs
#Affy
#data_norm_1 = read.table("geneset_1.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_2 = read.table("geneset_2.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_3 = read.table("geneset_3.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand_1 = read.table("geneset_1_2_3_random_2.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#Affy non-zero
data_norm_1 = read.table("geneset_1.Affy_by_clustersize_nonzero.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_2 = read.table("geneset_2.Affy_by_clustersize_nonzero.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_3 = read.table("geneset_3.Affy_by_clustersize_nonzero.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_rand_1 = read.table("geneset_1_2_3_random_2.Affy_by_clustersize_nonzero.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)


#Specify a maximum cluster size
maxclustersize=2
#maxclustersize=10

data_norm_1_undermax = data_norm_1[data_norm_1[,1]<=maxclustersize,][,2]
data_norm_2_undermax = data_norm_2[data_norm_2[,1]<=maxclustersize,][,2]
data_norm_3_undermax = data_norm_3[data_norm_3[,1]<=maxclustersize,][,2]
data_rand_1_undermax = data_rand_1[data_rand_1[,1]<=maxclustersize,][,2]

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
results_1=cumfreq(data_norm_1_undermax)
results_2=cumfreq(data_norm_2_undermax)
results_3=cumfreq(data_norm_3_undermax)
results_4=cumfreq(data_rand_1_undermax)


#Prepare to plot cum fractions for each on one plot
#plot(c(results_1[,"scores"],results_2[,"scores"],results_3[,"scores"]), c(results_1[,"cum_freqs"],results_2[,"cum_freqs"],results_3[,"cum_freqs"]), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Coexpression score (Pearson)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
#plot(c(results_1[,"scores"],results_2[,"scores"],results_3[,"scores"]), c(results_1[,"cum_freqs"],results_2[,"cum_freqs"],results_3[,"cum_freqs"]), type="n", xlim=c(0,1), ylim=c(0.75,1), xlab="Coexpression score (Pearson)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
#plot(c(results_1[,"scores"],results_2[,"scores"],results_3[,"scores"],results_4[,"scores"]), c(results_1[,"cum_freqs"],results_2[,"cum_freqs"],results_3[,"cum_freqs"],results_4[,"cum_freqs"]), type="n", xlim=c(0,1), ylim=c(0.75,1), xlab="Coexpression score (Pearson)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
plot(c(results_1[,"scores"],results_2[,"scores"],results_3[,"scores"],results_4[,"scores"]), c(results_1[,"cum_freqs"],results_2[,"cum_freqs"],results_3[,"cum_freqs"],results_4[,"cum_freqs"]), type="n", xlim=c(-1,1), ylim=c(0,1), xlab="Coexpression score (Pearson)", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)

#Add points
points(results_1[,"scores"], results_1[,"cum_freqs"], pch=17, cex=0.25, col="orange")
points(results_2[,"scores"], results_2[,"cum_freqs"], pch=19, cex=0.25, col="darkgreen")
points(results_3[,"scores"], results_3[,"cum_freqs"], pch=22, cex=0.25, col="blue")
points(results_4[,"scores"], results_4[,"cum_freqs"], pch=23, cex=0.25, col="red")

#Plot a legend
legend=c("genelist1", "genelist2", "genelist3", "random")
legend(0.4,0.4, legend, pch=c(17,19,22,23), col=c("orange","darkgreen","blue", "red"), cex=1.25)

dev.off()

#Test for significance between distributions
ks.test(x=coexp_norm_1_sorted, y=coexp_norm_2_sorted, alternative = c("two.sided"), exact = NULL)
ks.test(x=coexp_norm_1_sorted, y=coexp_norm_3_sorted, alternative = c("two.sided"), exact = NULL)
ks.test(x=coexp_norm_2_sorted, y=coexp_norm_3_sorted, alternative = c("two.sided"), exact = NULL)

