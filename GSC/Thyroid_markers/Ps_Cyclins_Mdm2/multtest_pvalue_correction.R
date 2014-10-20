############################################
library(multtest)
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Ps_Cyclins_Mdm2")
#pvalue_data=read.table("marker_vs_diag_pvalues.txt",header=TRUE,na.strings="NA",sep="\t")
pvalue_data=read.table("marker_vs_prog_pvalues.txt",header=TRUE,na.strings="NA",sep="\t")


Cont_pvalues=pvalue_data

cont_pvalues_adj=mt.rawp2adjp(Cont_pvalues[,3], proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY"))

cont_pvalues_adj_orig_order = cont_pvalues_adj$adjp[order(cont_pvalues_adj$index),]

#write.table (cont_pvalues_adj_orig_order, file="marker_vs_diag_pvalues_adj.txt", quote=FALSE, sep="\t")
write.table (cont_pvalues_adj_orig_order, file="marker_vs_prog_pvalues_adj.txt", quote=FALSE, sep="\t")
