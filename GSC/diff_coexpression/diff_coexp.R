datadir = "/home/obig/Projects/diff_coexpression/R_scripts"
dir(datadir)
setwd(datadir)

#file1 = read.table("sage_common_normal_zeros.matrix", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
#file2 = read.table("sage_common_cancer_zeros.matrix", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

file1 = read.table("sage_common_normal_nulls.matrix", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
file2 = read.table("sage_common_cancer_nulls.matrix", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#file1 = read.table("file1.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
#file2 = read.table("file2.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

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
write (file="corr_results.2.txt", t(output), ncolumns=2)

