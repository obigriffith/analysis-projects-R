datadir = "~/Projects/cisRED_coexpression/MotifDistanceAnalysis/motif_distance_files"
setwd(datadir)

data_coexp = read.table("coexpressed_TMMgt7_allgenes.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_coexp = read.table("coexpressed_multiplatform.out2", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_noncoexp = read.table("random_TMMgt7_allgenes.out", header=F, quote="", sep="\t", comment.char="", as.is=1)

resultdir="~/Projects/cisRED_coexpression/MotifDistanceAnalysis/R_analysis"
setwd(resultdir)
#Postscript file
ps_file="Motifdist_TMM_allgenes_coexp_vs_random.ps"
#ps_file="Motifdist_multiplatform_coexp_vs_noncoexp.ps"

#Filter down to only entries where motif distance is valid (distance of -1 is used to indicated one or more gene in the pair did not have motifs in the version of cisred used)
data_coexp_valid=data_coexp[data_coexp[,4]>=0,]
data_noncoexp_valid=data_noncoexp[data_noncoexp[,4]>=0,]

coexp_motif_dists=data_coexp_valid[,4]
noncoexp_motif_dists=data_noncoexp_valid[,4]

hist_coexp = hist(coexp_motif_dists,plot=FALSE)
hist_noncoexp = hist(noncoexp_motif_dists,plot=FALSE)


#Calculate t-test between coexpressed an non-coexpressed
coexp_vs_noncoexp_ttest=t.test(x=coexp_motif_dists,y=noncoexp_motif_dists, alternative=c("two.sided"), mu=0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
coexp_vs_noncoexp_ttest_pvalue=coexp_vs_noncoexp_ttest$p.value
coexp_vs_noncoexp_ttest_meanx = coexp_vs_noncoexp_ttest$estimate[1]
coexp_vs_noncoexp_ttest_meany = coexp_vs_noncoexp_ttest$estimate[2]


#Plot coexpressed motif distances versus non-coexpressed distances
postscript(ps_file, pointsize=1, width=9.5, height=7)
plot(c(hist_coexp$mids,hist_noncoexp$mids), c(hist_coexp$counts,hist_noncoexp$counts), xlab="Minimum Motif distance", ylab="counts", main="distribution of minimum motif distance for KSL coexpressed genes vs non-coexpressed",cex.lab=2.5, cex.axis=2.5, cex.main=2)
points(hist_coexp$mids, hist_coexp$counts, pch=19, cex=2.5)
points(hist_noncoexp$mids, hist_noncoexp$counts, pch=23, cex=2.5)
legend=c("Coexpressed","Non-Coexpressed")
legend(0.6,60, legend, pch=c(19,23),cex=2.5)

#text(0.6, 40, c("pvalue=",coexp_vs_noncoexp_ttest_pvalue,"mean_coexp=",coexp_vs_noncoexp_ttest_meanx,"mean_noncoexp=",coexp_vs_noncoexp_ttest_meany))

dev.off()