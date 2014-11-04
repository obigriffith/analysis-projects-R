datadir = "/home/obig/Projects/diff_coexpression/data/prostate/Singh_prostate/"

dir(datadir)
setwd(datadir)

file1 = read.table("Prostate_TN_final0701_allmeanScale.Normal", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
file2 = read.table("Prostate_TN_final0701_allmeanScale.Tumor", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#Switch files for tumor vs normal (instead of normal vs tumor).
#file2 = read.table("Prostate_TN_final0701_allmeanScale.Normal", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
#file1 = read.table("Prostate_TN_final0701_allmeanScale.Tumor", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

n=10 #Set this as the neighborhood to examine.

#Find all pairwise correlations for each file
#file1_cor = cor(t(file1), method = "pearson", use = "pairwise.complete.obs")
#file2_cor = cor(t(file2), method = "pearson", use = "pairwise.complete.obs")
file1_cor = cor(t(file1), method = "spearman", use = "pairwise.complete.obs")
file2_cor = cor(t(file2), method = "spearman", use = "pairwise.complete.obs")

num_genes = length(file1_cor[,1]) #Determines number of genes to loop through

#Instead try a simple loop
for (i in 1:num_genes){
file1_cor[i,i]=NA
file2_cor[i,i]=NA
}

#Create vector to store genenames
gene_names = dimnames(file1)[[1]]

#Create a matrix to store correlation of correlations between genes
#gene_cors_cors = array(0, dimnames = list(gene_names, "cor_cor", dim=c(length(gene_names),1))
gene_cors_cors = array(0, dimnames = list(gene_names, c("cor_cor", "normal_mean","tumor_mean","normal_var","tumor_var","mean_diff","ttest_pval")), dim=c(length(gene_names),7))

#For a specified gene get closest n neighbors from normal and corresponding distances in tumor
#my_gene='928_at'
#my_gene='41706_at'
#my_gene='41661_at'
#my_gene='2041_i_at'
#my_gene='41480_at'
#my_gene='1874_at'
#my_gene='33436_at'
my_gene='2038_g_at'

norm=file1_cor[,my_gene]   #vector of coexpression distances for gene i in normal data
tumor=file2_cor[,my_gene]   #vector of coexpression distances for gene i in tumour data
norm_ordered=order(norm,decreasing=TRUE)   #order returns a vector of indices for genes sorted by value.  ie the first entry is the position in the unsorted vector for the highest value
norm_top_indices=norm_ordered[1:n]   #indices of top n most similar genes in normal vector
norm_ranked=rank(-norm,na.last=TRUE)  #rank returns vector of ranks for each gene pair.  Note default order is increasing (no decreasing option).  For numerical vector we can multiply all values by -1 to get around this
tumor_ranked=rank(-tumor,na.last=TRUE) #do the same for tumor
top_ranks_in_norm=norm_ranked[norm_top_indices] #Get ranks for top genes as determined by the norm_ordered vector
top_ranks_in_tumor=tumor_ranked[norm_top_indices] #do the same for tumor
norm_top_dist=norm[norm_top_indices]   #distances of top n most similar genes to gene i in normal vector
tumor_top_dist=tumor[norm_top_indices]   #distances of corresponding genes to gene i in tumor vector

norm_top_genes=labels(norm_top_dist) #Get labels (ie. probe_ids) for top n genes to my_gene.  These should be identical for norm_top_dist and tumor_top_dist
my_gene_norm_top_genes=c(my_gene,norm_top_genes)

#Get expression values for gene of interest and top n neighbors in norm and tumor
norm_my_gene_top_genes_exp=file1[my_gene_norm_top_genes,]
tumor_my_gene_top_genes_exp=file2[my_gene_norm_top_genes,]


resultdir = "/home/obig/Projects/diff_coexpression/results/prostate/Singh_prostate/Spearman/top_gene_profiles/"
setwd(resultdir)
output = cbind(my_gene_norm_top_genes,norm_my_gene_top_genes_exp)
write (file="Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_top10_norm_profiles_2038_g_at.txt", t(output), ncolumns=51)
output = cbind(my_gene_norm_top_genes,tumor_my_gene_top_genes_exp)
write (file="Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_top10_tumor_profiles_2038_g_at.txt", t(output), ncolumns=53)

#For tumor vs normal
#output = cbind(my_gene_norm_top_genes,norm_my_gene_top_genes_exp)
#write (file="Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_top10_tumor_profiles_41480_at.txt", t(output), ncolumns=53)
#output = cbind(my_gene_norm_top_genes,tumor_my_gene_top_genes_exp)
#write (file="Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_top10_norm_profiles_41480_at.txt", t(output), ncolumns=51)




