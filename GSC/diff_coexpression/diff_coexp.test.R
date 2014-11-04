datadir = "/home/obig/Projects/diff_coexpression/R_scripts"
dir(datadir)
setwd(datadir)

n=5 #Set this as the neighborhood to examine.

file1 = read.table("test1.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
file2 = read.table("test2.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#Find all pairwise correlations for each file
#file1_cor = cor(t(file1), method = "pearson", use = "pairwise.complete.obs")
#file2_cor = cor(t(file2), method = "pearson", use = "pairwise.complete.obs")
file1_cor = cor(t(file1), method = "pearson", use = "complete.obs")
file2_cor = cor(t(file2), method = "pearson", use = "complete.obs")

num_genes = length(file1_cor[,1]) 

#Remove self-correlations (ie. geneA versus geneA, etc.)
diag(file1_cor) = NA
diag(file2_cor) = NA

#Instead try a simple loop
#for (i in 1:num_genes){
#file1_cor[i,i]=0
#file2_cor[i,i]=0
#}

#Create vector to store genenames
gene_names = dimnames(file1)[[1]]

#Create a matrix to store correlation of correlations between genes
#gene_cors_cors = array(0, dimnames = list(gene_names, "cor_cor"), dim=c(length(gene_names),1))
gene_cors_cors = array(0, dimnames = list(gene_names, c("cor_cor", "normal_mean","tumor_mean","norm_var","tumor_var","mean_diff","ttest_pval")), dim=c(length(gene_names),7))


#For each gene get top n neighbors from normal and corresponding distances in tumor
for (i in 1:num_genes){
   norm=file1_cor[,i]   #vector of coexpression distances for gene i in normal data
   tumor=file2_cor[,i]   #vector of coexpression distances for gene i in tumour data
   norm_ordered=order(norm,decreasing=TRUE)   #order returns a vector of indices for genes sorted by value.  ie the first entry is the position in the unsorted vector for the lowest value
   norm_top_indices=norm_ordered[1:n]   #indices of top n most similar genes in normal vector
   norm_ranked=rank(-norm,na.last=TRUE)  #rank returns
   tumor_ranked=rank(-tumor,na.last=TRUE)
   top_ranks_in_norm=norm_ranked[norm_top_indices]
   top_ranks_in_tumor=tumor_ranked[norm_top_indices]
   norm_top_dist=norm[norm_top_indices]   #distances of top n most similar genes to gene i in normal vector
   tumor_top_dist=tumor[norm_top_indices]   #distances of corresponding genes to gene i in tumor vector
   #Calculate statistics about gene
   cor_cor = cor(x=top_ranks_in_norm, y=top_ranks_in_tumor, method = "pearson", use = "complete.obs")
   norm_top_dist_mean = mean(norm_top_dist, trim = 0, na.rm = TRUE)
   tumor_top_dist_mean = mean(tumor_top_dist, trim = 0, na.rm = TRUE)
   norm_top_dist_var = var(norm_top_dist, y = NULL, na.rm = TRUE)
   tumor_top_dist_var = var(tumor_top_dist, y = NULL, na.rm = TRUE)
   norm_vs_tumor_ttest=t.test(x=norm_top_dist,y=tumor_top_dist,alternative = c("less"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
   norm_vs_tumor_meandiff=norm_top_dist_mean - tumor_top_dist_mean

#   cor_cor = cor(x=norm_top_dist, y=tumor_top_dist, method = "pearson", use = "complete.obs")
#   cor_cor = cor(x=norm_top_dist, y=tumor_top_dist, method = "pearson", use = "pairwise.complete.obs")
   gene_cors_cors[i,"cor_cor"] = cor_cor
   gene_cors_cors[i,"normal_mean"] = norm_top_dist_mean
   gene_cors_cors[i,"tumor_mean"] = tumor_top_dist_mean
   gene_cors_cors[i,"norm_var"] = norm_top_dist_var
   gene_cors_cors[i,"tumor_var"] = tumor_top_dist_var
   gene_cors_cors[i,"mean_diff"] = norm_vs_tumor_meandiff
   gene_cors_cors[i,"ttest_pval"] = norm_vs_tumor_ttest$p.value
}

output = cbind(gene_names,gene_cors_cors)
write (file="test_output.txt", t(output), ncolumns=8)

#Produce graph of distribution
postscript("test_hist.2.ps", pointsize=1)
#hist(gene_cors_cors, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="Correlation of Pearson Correlations for all genes comparing normal and tumor", plot=TRUE)
#hist(gene_cors_cors[,"cor_cor"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="Correlation of Pearson Correlations for all genes comparing normal and tumor", plot=TRUE)
