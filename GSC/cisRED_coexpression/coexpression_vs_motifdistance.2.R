#datadir = "~/Projects/cisRED_coexpression/MotifDistanceAnalysis/motif_distance_files"
#datadir = "~/Projects/cisRED_coexpression/MotifDistanceAnalysis/motif_distance_files/cisred_2/pcutoff_25"
datadir = "~/Projects/cisRED_coexpression/MotifDistanceAnalysis/motif_distance_files/cisred_2/pcutoff_05"
setwd(datadir)

#cisred_1_2e
#data_coexp = read.table("coexpressed_multiplatform_allgenes.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_coexp = read.table("coexpressed_TMMgt7_allgenes.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_coexp = read.table("coexpressed_TMMgt7.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_coexp = read.table("coexpressed_multiplatform.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_noncoexp = read.table("noncoexpressed_top10.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_noncoexp = read.table("random_TMMgt7_allgenes.out", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_noncoexp = read.table("random_multiplatform_allgenes.out", header=F, quote="", sep="\t", comment.char="", as.is=1)

#cisred_2 - all genes
data_coexp = read.table("coexpressed_multiplatform.clean", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_coexp = read.table("coexpressed_TMMgt7.clean", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_noncoexp = read.table("random_multiplatform.out", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_noncoexp = read.table("random_TMMgt7.out", header=F, quote="", sep="\t", comment.char="", as.is=1)

resultdir="~/Projects/cisRED_coexpression/MotifDistanceAnalysis/R_analysis/cisred_2/pcutoff_05"
setwd(resultdir)
#Postscript file
#ps_file="Motifdist_TMM_coexp_vs_noncoexp.ps"
#ps_file="Motifdist_multiplatform_coexp_vs_noncoexp.ps"
#ps_file="Motifdist_TMM_all_genes_coexp_vs_noncoexp.ps"
#ps_file="Motifdist_multiplatform_all_genes_coexp_vs_noncoexp.ps"
#ps_file="Motifdist_multiplatform_all_genes_coexp_vs_random_cumul.ps"
#ps_file="Motifdist_TMM_all_genes_coexp_vs_noncoexp_cumul.ps"
pdf_file="Motifdist_multiplatform_all_genes_coexp_vs_random_cumul.pdf"
#pdf_file="Motifdist_TMM_all_genes_coexp_vs_random_cumul.pdf"


#Filter down to only entries where motif distance is valid (distance of -1 is used to indicated one or more gene in the pair did not have motifs in the version of cisred used)
data_coexp_valid=data_coexp[data_coexp[,4]>=0,]
data_noncoexp_valid=data_noncoexp[data_noncoexp[,4]>=0,]

coexp_motif_dists=data_coexp_valid[,4]
noncoexp_motif_dists=data_noncoexp_valid[,4]

#Calculate t-test between coexpressed and non-coexpressed
coexp_vs_noncoexp_ttest=t.test(x=coexp_motif_dists,y=noncoexp_motif_dists, alternative=c("two.sided"), mu=0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
coexp_vs_noncoexp_ttest_pvalue=coexp_vs_noncoexp_ttest$p.value
coexp_vs_noncoexp_ttest_meanx = coexp_vs_noncoexp_ttest$estimate[1]
coexp_vs_noncoexp_ttest_meany = coexp_vs_noncoexp_ttest$estimate[2]

#Create histogram without plots for later plotting on the same graph
hist_coexp = hist(coexp_motif_dists,plot=FALSE, breaks=20)
hist_noncoexp = hist(noncoexp_motif_dists,plot=FALSE, breaks=20)

#calculate fractions and cumulative fractions for each bin
#hist_coexp_fractions = hist_coexp$counts/sum(hist_coexp$counts)
#hist_noncoexp_fractions = hist_noncoexp$counts/sum(hist_noncoexp$counts)
hist_coexp_cum_fractions=cumsum(hist_coexp$counts/sum(hist_coexp$counts))
hist_noncoexp_cum_fractions=cumsum(hist_noncoexp$counts/sum(hist_noncoexp$counts))

#Plot coexpressed motif distances versus non-coexpressed distances
#postscript(ps_file, pointsize=1, width=9.0, height=7)
pdf(pdf_file, width=10.5, height=8)

plot(c(hist_coexp$mids,hist_noncoexp$mids), c(hist_coexp_cum_fractions,hist_noncoexp_cum_fractions), type="n", xlim=range(0,1), xlab="Minimum Motif Pair Distance", ylab="Cumulative Fraction of Motif Pair Distances", main="Distribution of minimum motif distance for coexpressed genes versus random genes",cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
points(hist_coexp$mids, hist_coexp_cum_fractions, pch=19,col=1, cex=1)
points(hist_noncoexp$mids, hist_noncoexp_cum_fractions, pch=22,col=2, cex=1)
#lines(lowess(x=hist_coexp$mids, y=hist_coexp_cum_fractions), lty=1,col = 1)
#lines(lowess(x=hist_noncoexp$mids, y=hist_noncoexp_cum_fractions), lty=1, col = 2)
lines(x=hist_coexp$mids, y=hist_coexp_cum_fractions, lty=1,col = 1)
lines(x=hist_noncoexp$mids, y=hist_noncoexp_cum_fractions, lty=1, col = 2)
legend=c("Coexpressed","Random")
legend(0.7,0.2, legend, lty=c(1,1), pch=c(19,22),col = c(1,2),cex=1.4)
#legend(0.7,0.1, legend, lty=c(1,1), pch=c(19,22),col = c(1,2),cex=2.5)
#text(0.6, 40, c("pvalue=",coexp_vs_noncoexp_ttest_pvalue,"mean_coexp=",coexp_vs_noncoexp_ttest_meanx,"mean_noncoexp=",coexp_vs_noncoexp_ttest_meany))

dev.off()
