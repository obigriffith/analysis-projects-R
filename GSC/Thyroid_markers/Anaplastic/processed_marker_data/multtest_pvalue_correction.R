library(multtest)
pvalue_data_group1=read.table("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_group1.txt",header=TRUE,sep="\t",na.strings="NA")
pvalue_data_group2=read.table("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_group2.txt",header=TRUE,sep="\t",na.strings="NA")
pvalue_data_MH=read.table("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_MH.txt",header=TRUE,sep="\t",na.strings="NA")


Cont_pvalues_grouping1=pvalue_data_group1[,2]
Cont_pvalues_grouping2=pvalue_data_group2[,2]
MH_pvalues=pvalue_data_MH[,2]

cont_pvalues_grouping1_adj=mt.rawp2adjp(Cont_pvalues_grouping1, proc=c("Bonferroni","BH"))
cont_pvalues_grouping2_adj=mt.rawp2adjp(Cont_pvalues_grouping2, proc=c("Bonferroni","BH"))
MH_pvalues_adj=mt.rawp2adjp(MH_pvalues, proc=c("Bonferroni","BH"))

cont_pvalues_grouping1_adj_orig_order = cont_pvalues_grouping1_adj$adjp[order(cont_pvalues_grouping1_adj$index),]
cont_pvalues_grouping2_adj_orig_order = cont_pvalues_grouping2_adj$adjp[order(cont_pvalues_grouping2_adj$index),]
MH_pvalues_adj_orig_order = MH_pvalues_adj$adjp[order(MH_pvalues_adj$index),]

#result_matrix=cbind(cont_pvalues_grouping1_adj_orig_order,cont_pvalues_grouping2_adj_orig_order,MH_pvalues_adj_orig_order)
#row.names(result_matrix)=pvalue_data[,1]

row.names(cont_pvalues_grouping1_adj_orig_order)=pvalue_data_group1[,1]
row.names(cont_pvalues_grouping2_adj_orig_order)=pvalue_data_group2[,1]
row.names(MH_pvalues_adj_orig_order)=pvalue_data_MH[,1]

write.table (cont_pvalues_grouping1_adj_orig_order, sep="\t", file="C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_group1_adj.txt")
write.table (cont_pvalues_grouping2_adj_orig_order, sep="\t", file="C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_group2_adj.txt")
write.table (MH_pvalues_adj_orig_order, sep="\t", file="C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/pvalues_MH_adj.txt")
