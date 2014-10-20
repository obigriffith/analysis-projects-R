library(multtest)

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide"
setwd(outdir)

#genes of interest
gene_interest_file="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide/HumanGlycosyltransferases_etc.txt"
gene_interest_import=read.table(gene_interest_file, header = TRUE, na.strings = "NA", sep="\t", as.is=1:2)
genes_interest=gene_interest_import[,1]
rownames(gene_interest_import)=gene_interest_import[,1]

#Cell line info
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide/BCCL_data_list.txt"
cell_data_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=1)
lines=cell_data_import[,1]

### exonarray data ###
EAdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ExonArray/prefiltered/breastExon_genelevel.csv"
EA_data_import=read.csv(EAdatafile) 
EA_feat_data=as.vector(EA_data_import[,1])
EA_data=EA_data_import[,2:length(colnames(EA_data_import))]

#Get only genes and cell lines of interest
EA_data=EA_data[EA_feat_data %in% genes_interest,lines]
EA_feat_data=EA_feat_data[EA_feat_data %in% genes_interest]

#Set rownames for EA data to gene symbol. Should be unique now
rownames(EA_data)=EA_feat_data

#Retrieve corresponding list name
list_name=gene_interest_import[EA_feat_data,2]

#Calculate means for TN (first 5 lines) versus Luminal (last 5 lines) and FC
TN_means=apply(EA_data[,1:5],1,mean)
Lum_means=apply(EA_data[,6:10],1,mean)
FC=TN_means/Lum_means
FC[which(FC<1)]= -1/FC[which(FC<1)]


#Unlog2 data and calculate FCs
EA_data_raw=2^EA_data
TN_means_raw=apply(EA_data_raw[,1:5],1,mean)
Lum_means_raw=apply(EA_data_raw[,6:10],1,mean)
FC_raw=TN_means_raw/Lum_means_raw
FC_raw[which(FC_raw<1)]= -1/FC_raw[which(FC_raw<1)]

#calculate a wilcox p-value
wilcox_results=rep(NA,length(rownames(EA_data)))
for (i in 1:length(rownames(EA_data))){
  p=wilcox.test(x=as.numeric(EA_data[i,1:5]), y=as.numeric(EA_data[i,6:10]), alternative="two.sided")$p.value
  wilcox_results[i]=p
}

#Correct p-values
pvalues=as.numeric(wilcox_results)
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]

#Write results to file
EA_outdata=cbind(rownames(EA_data),list_name,TN_means,Lum_means,FC,FC_raw,wilcox_results,pvalues_adj_orig_order[,2],EA_data)
write.table(EA_outdata, file="GlycoProteins_exonarray_analysis.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

