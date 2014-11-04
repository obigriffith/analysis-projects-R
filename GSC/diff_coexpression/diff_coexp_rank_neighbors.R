datadir = "/home/obig/Projects/diff_coexpression/data/prostate/Singh_prostate/"
dir(datadir)
setwd(datadir)

n=5 #Set this as the neighborhood to examine.

file1 = read.table("Prostate_TN_final0701_allmeanScale.Normal", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
file2 = read.table("Prostate_TN_final0701_allmeanScale.Tumor", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)



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

#Create vector to store genenames
gene_names = dimnames(file1)[[1]]

#Create a matrix to store correlation of correlations between genes
gene_cors_cors = array(0, dimnames = list(gene_names, "cor_cor"), dim=c(length(gene_names),1))

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
   cor_cor = cor(x=top_ranks_in_norm, y=top_ranks_in_tumor, method = "pearson", use = "complete.obs") #Finally, calculate a Pearson correlation of the ranks of pearson correlations (same as Spearman of Pearsons)
   gene_cors_cors[i,"cor_cor"] = cor_cor  #Enter result into matrix for printout
}

output = cbind(gene_names,gene_cors_cors)
write (file="Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top5_diff_coexp_byrank.txt", t(output), ncolumns=2)

#Produce graph of distribution
postscript("Prostate_TN_final0701_allmeanScale.Normal_vs_Tumor_top5_byrank_freq_dist.ps", pointsize=1)
hist(gene_cors_cors, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="Spearman Correlation of Pearson Correlations for all genes comparing normal (closest 10 neighbors) and tumor", plot=TRUE)

