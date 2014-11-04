datadir = "/home/obig/Projects/diff_coexpression/data/prostate/Singh_prostate/"
dir(datadir)
setwd(datadir)

file1 = read.table("Prostate_TN_final0701_allmeanScale.Normal", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
file2 = read.table("Prostate_TN_final0701_allmeanScale.Tumor", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#Find all pairwise correlations for each file
file1_cor = cor(t(file1), method = "pearson", use = "pairwise.complete.obs")
file2_cor = cor(t(file2), method = "pearson", use = "pairwise.complete.obs")


gene_names = dimnames(file1)[[1]]
gene_pear_cors = array(0, dimnames = list(gene_names, "pearson_cor"), dim=c(length(gene_names),1))

num_genes = length(file1_cor[,1]) 


for (i in 1:num_genes){

   gene_cor_cor = cor(x=file1_cor[,i], y=file2_cor[,i], method = "pearson", use = "pairwise.complete.obs")
   gene_pear_cors[i,"pearson_cor"] = gene_cor_cor
   #print(gene_cor_cor)
}

output = cbind(gene_names,gene_pear_cors)
write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_diff_coexp.txt", t(output), ncolumns=2)

#Produce graph of distribution
postscript("Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_freq_dist.ps", pointsize=1)
hist(gene_pear_cors, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="Correlation of Pearson Correlations for all genes comparing normal and tumor", plot=TRUE)
