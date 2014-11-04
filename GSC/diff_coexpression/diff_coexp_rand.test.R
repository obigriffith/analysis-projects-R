datadir = "/home/obig/Projects/diff_coexpression/R_scripts"
dir(datadir)
setwd(datadir)

n=10 #Set this as the neighborhood to examine.

file1 = read.table("test1.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)
#file2 = read.table("test2.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#For random permutations we want to create two pseudo files from either the normal set or tumor set
#First determine the number of columns and put in a vector
exp_count=length(file1[1,])
exp_indices=(1:exp_count)

#now randomly sample from the vector n times with replacement (where n=number of experiments)
#Create two permutations of the indices
perm_ind1=sample(exp_indices, exp_count, replace = TRUE, prob = NULL)
perm_ind2=sample(exp_indices, exp_count, replace = TRUE, prob = NULL)

file1_perm1=file1[perm_ind1]
file1_perm2=file1[perm_ind2]

#Find all pairwise correlations for each file
file1_cor = cor(t(file1_perm1), method = "pearson", use = "pairwise.complete.obs")
file2_cor = cor(t(file1_perm2), method = "pearson", use = "pairwise.complete.obs")

#Remove self-correlations (ie. geneA versus geneA, etc.)
diag(file1_cor) = NA
diag(file2_cor) = NA

#Create vector to store genenames
gene_names = dimnames(file1)[[1]]

#Create a matrix to store correlation of correlations between genes
gene_cors_cors = array(0, dimnames = list(gene_names, "cor_cor"), dim=c(length(gene_names),1))

#For each gene get top n neighbors from normal and corresponding distances in tumor
num_genes = length(file1_cor[,1]) 
for (i in 1:num_genes){
   norm=file1_cor[,i]   #vector of coexpression distances for gene i in normal data
   tumor=file2_cor[,i]   #vector of coexpression distances for gene i in tumour data
   norm_ordered=order(norm,decreasing=TRUE)   #order returns a vector of indices for genes sorted by value.  ie the first entry is the position in the unsorted vector for the highest value
   norm_top_indices=norm_ordered[1:n]   #indices of top n most similar genes in normal vector
   norm_top_dist=norm[norm_top_indices]   #distances of top n most similar genes to gene i in normal vector
   tumor_top_dist=tumor[norm_top_indices]   #distances of corresponding genes to gene i in tumor vector
   cor_cor = cor(x=norm_top_dist, y=tumor_top_dist, method = "pearson", use = "pairwise.complete.obs")
   gene_cors_cors[i,"cor_cor"] = cor_cor
}

output = cbind(gene_names,gene_cors_cors)
write (file="test_output.txt", t(output), ncolumns=2)

#Produce graph of distribution
postscript("test_hist.ps", pointsize=1)
hist(gene_cors_cors, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="Correlation of Pearson Correlations for all genes comparing normal and tumor", plot=TRUE)
