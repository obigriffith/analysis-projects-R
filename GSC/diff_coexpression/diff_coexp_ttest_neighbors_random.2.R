datadir = "/home/obig/Projects/diff_coexpression/data/prostate/Singh_prostate/"
dir(datadir)
setwd(datadir)

n=10 #Set this as the neighborhood to examine.

infile1 = read.table("Prostate_TN_final0701_allmeanScale.Normal", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
#infile1 = read.table("Prostate_TN_final0701_allmeanScale.Tumor", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#For randomizations.  Just read in one file and then make two pseudofiles from it by sampling the columns with replacement
num_cols=length(infile1[1,])
file1 = sample(infile1, num_cols, replace=TRUE)
file2 = sample(infile1, num_cols, replace=TRUE)

#Find all pairwise correlations for each file
file1_cor = cor(t(file1), method = "pearson", use = "pairwise.complete.obs")
file2_cor = cor(t(file2), method = "pearson", use = "pairwise.complete.obs")

num_genes = length(file1_cor[,1]) #Determines number of genes to loop through

#Remove self-correlations (ie. geneA versus geneA, etc.)
#For some stupid reason, diag is giving memory errors on large matrices
#diag(file1_cor) = NA
#diag(file2_cor) = NA

#Instead try a simple loop
for (i in 1:num_genes){
file1_cor[i,i]=NA
file2_cor[i,i]=NA
}

#Would be good to graph distribution of Pearson's for each file - might be memory problems?
#postscript("Prostate_TN_final0701_allmeanScale.Tumor_Pearsons_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"ttest_pval"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="P-value", main="Tumor Pearson values for all genes pairs", plot=TRUE)
#postscript("Prostate_TN_final0701_allmeanScale.Normal_Pearsons_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"ttest_pval"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="P-value", main="Normal Pearson values for all genes pairs", plot=TRUE)

#Create vector to store genenames
gene_names = dimnames(file1)[[1]]

#Create a matrix to store correlation of correlations between genes
#gene_cors_cors = array(0, dimnames = list(gene_names, "cor_cor", dim=c(length(gene_names),1))
#gene_cors_cors = array(0, dimnames = list(gene_names, c("cor_cor", "normal_mean","tumor_mean","normal_var","tumor_var","mean_diff","ttest_pval")), dim=c(length(gene_names),7))
gene_cors_cors = array(0, dimnames = list(gene_names, "mean_diff"), dim=c(length(gene_names),1))

#For each gene get top n neighbors from normal and corresponding distances in tumor
for (i in 1:num_genes){
   norm=file1_cor[,i]   #vector of coexpression distances for gene i in normal data
   tumor=file2_cor[,i]   #vector of coexpression distances for gene i in tumour data
   norm_ordered=order(norm,decreasing=TRUE)   #order returns a vector of indices for genes sorted by value.  ie the first entry is the position in the unsorted vector for the highest value
   norm_top_indices=norm_ordered[1:n]   #indices of top n most similar genes in normal vector
   norm_ranked=rank(-norm,na.last=TRUE)  #rank returns vector of ranks for each gene pair.  Note default order is increasing (no decreasing option).  For numerical vector we can multiply all values by -1 to get around this
   tumor_ranked=rank(-tumor,na.last=TRUE) #do the same for tumor
   top_ranks_in_norm=norm_ranked[norm_top_indices] #Get ranks for top genes as determined by the norm_ordered vector
   top_ranks_in_tumor=tumor_ranked[norm_top_indices] #do the same for tumor
   norm_top_dist=norm[norm_top_indices]   #distances of top n most similar genes to gene i in normal vector
   tumor_top_dist=tumor[norm_top_indices]   #distances of corresponding genes to gene i in tumor vector
   #Calculate statistics about gene
#   cor_cor = cor(x=top_ranks_in_norm, y=top_ranks_in_tumor, method = "pearson", use = "complete.obs") #Finally, calculate a Pearson correlation of the ranks of pearson correlations (same as Spearman of Pearsons)
   norm_top_dist_mean = mean(norm_top_dist, trim = 0, na.rm = TRUE)
   tumor_top_dist_mean = mean(tumor_top_dist, trim = 0, na.rm = TRUE)
#   norm_top_dist_var = var(norm_top_dist, y = NULL, na.rm = TRUE)
#   tumor_top_dist_var = var(tumor_top_dist, y = NULL, na.rm = TRUE)
   norm_vs_tumor_meandiff=norm_top_dist_mean - tumor_top_dist_mean
#   norm_vs_tumor_ttest=t.test(x=norm_top_dist,y=tumor_top_dist,alternative = c("greater"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
   #fill summary list with all values
#   gene_cors_cors[i,"cor_cor"] = cor_cor  #Enter result into matrix for printout
#   gene_cors_cors[i,"normal_mean"] = norm_top_dist_mean
#   gene_cors_cors[i,"tumor_mean"] = tumor_top_dist_mean
#   gene_cors_cors[i,"normal_var"] = norm_top_dist_var
#   gene_cors_cors[i,"tumor_var"] = tumor_top_dist_var
   gene_cors_cors[i,"mean_diff"] = norm_vs_tumor_meandiff
#   gene_cors_cors[i,"ttest_pval"] = norm_vs_tumor_ttest$p.value
}

resultdir = "/home/obig/Projects/diff_coexpression/results/prostate/Singh_prostate/top10/random/"
setwd(resultdir)

#Write summary results to file
#headers=as.matrix(c("probe","cor_cor", "normal_mean","tumor_mean","norm_var","tumor_var","mean_diff","ttest_pval"))
#write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top10_diff_coexp_byrank.txt", t(headers), ncolumns=8)
#output = cbind(gene_names,gene_cors_cors)
#write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top10_diff_coexp_byrank.txt", t(output), ncolumns=8, append=TRUE)

headers=as.matrix(c("probe","mean_diff"))
write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top10_diff_coexp_byrank.test.txt", t(headers), ncolumns=2)
output = cbind(gene_names,gene_cors_cors)
write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top10_diff_coexp_byrank.test.txt", t(output), ncolumns=2, append=TRUE)

#Produce graph of distribution of correlation of correlations
#postscript("Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top10_byrank_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"cor_cor"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Spearman correlation coefficient (r)", main="Spearman Correlation of Pearson Correlations for all genes comparing normal (closest 10 neighbors) and tumor", plot=TRUE)

#Produce graph of distribution of Normal Means
#postscript("Prostate_TN_final0701_allmeanScale.Normal_means_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"normal_mean"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Mean Pearson correlation coefficient (r)", main="Normal - Mean Pearson Correlations for all genes (closest 10 neighbors)", plot=TRUE)

#Produce graph of distribution of Tumor Means
#postscript("Prostate_TN_final0701_allmeanScale.Tumor_means_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"tumor_mean"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Mean Pearson correlation coefficient (r)", main="Tumor - Mean Pearson Correlations for all genes (closest 10 neighbors)", plot=TRUE)

#Produce graph of distribution of Normal Variances
#postscript("Prostate_TN_final0701_allmeanScale.Normal_vars_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"normal_var"], breaks=201,freq=TRUE, col="blue", xlim=range(-0.5,0.5), xlab="Variance of Pearson correlation coefficient (r)", main="Normal - Variance of Pearson Correlations for all genes (closest 10 neighbors)", plot=TRUE)

#Produce graph of distribution of Tumor Variances
#postscript("Prostate_TN_final0701_allmeanScale.Tumor_vars_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"tumor_var"], breaks=201,freq=TRUE, col="blue", xlim=range(-0.5,0.5), xlab="Variance of Pearson correlation coefficient (r)", main="Tumor - Mean Pearson Correlations for all genes (closest 10 neighbors)", plot=TRUE)

#Produce graph of distribution of Mean differences
#postscript("Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_mean_diffs_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"mean_diff"], breaks=201,freq=TRUE, col="blue", xlim=range(-1,2), xlab="Mean difference", main="Tumor vs Normal difference in Mean Pearson for all genes (closest 10 neighbors)", plot=TRUE)

#Produce graph of distribution of Pvalues
#postscript("Prostate_TN_final0701_allmeanScale.Tumor_vs_Normal_ttests_pvalue_freq_dist.ps", pointsize=1)
#hist(gene_cors_cors[,"ttest_pval"], breaks=201,freq=TRUE, col="blue", xlab="P-value", main="Tumor vs Normal t-test p-values for all genes (closest 10 neighbors)", plot=TRUE)
