library(multtest)
library(exactRankTests)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis")
datafile="BenignMalignant_58markers_23JAN08_LN_paired_ungrouped.txt"

outfile="BenignMalignant_58markers_23JAN08_ungrouped_LNanalysis.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

primary_marker_data=data[,36:92]
LN_marker_data=data[,93:149]

marker_names=colnames(primary_marker_data)

#Create array to store results for prognostic contingency table stats
LN_results = array(NA, dimnames = list(marker_names, c("Spearman_rho","Spearman_pvalue","Wilcox_exact_pvalue")), dim=c(length(marker_names),3))

for (i in 1:length(marker_names)){
 if(nlevels(as.factor(primary_marker_data[,i]))<2){
  next
 }
 if(nlevels(as.factor(LN_marker_data[,i]))<2){
  next
 }
 wilcox=wilcox.exact(x=primary_marker_data[,i], y=LN_marker_data[,i], alternative="two.sided", paired=TRUE)
 spearman=cor.test(x=primary_marker_data[,i], y=LN_marker_data[,i], use = "complete.obs", method = "spearman")
 LN_results[i,"Wilcox_exact_pvalue"]=wilcox$p.value
 LN_results[i,"Spearman_pvalue"]=spearman$p.value
 LN_results[i,"Spearman_rho"]=spearman$estimate
}

#Correct p-values
spearman_pvalues=as.numeric(LN_results[,"Spearman_pvalue"])
spearman_pvalues_adj=mt.rawp2adjp(spearman_pvalues, proc=c("Bonferroni","BH"))
spearman_pvalues_adj_orig_order=spearman_pvalues_adj$adjp[order(spearman_pvalues_adj$index),]
LN_results=cbind(LN_results, spearman_pvalues_adj_orig_order[,2:3])

wilcox_pvalues=as.numeric(LN_results[,"Wilcox_exact_pvalue"])
wilcox_pvalues_adj=mt.rawp2adjp(wilcox_pvalues, proc=c("Bonferroni","BH"))
wilcox_pvalues_adj_orig_order=wilcox_pvalues_adj$adjp[order(wilcox_pvalues_adj$index),]
LN_results=cbind(LN_results, wilcox_pvalues_adj_orig_order[,2:3])

#Write complete results to file
write.table (LN_results, sep="\t", file=outfile)
