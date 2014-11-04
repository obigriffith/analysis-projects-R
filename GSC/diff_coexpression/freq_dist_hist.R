datadir = "/projects/02/coexpression/distance_analysis/Gene_Matrix/human_gene_matrix/"
dir(datadir)
setwd(datadir)

file1 = read.table("human_affy_gene_matrix.txt", header=T, quote="", sep="\t", comment.char="", as.is=1, row.names=1)

#Find all pairwise correlations for each file
file1_cor = cor(t(file1), method = "pearson", use = "pairwise.complete.obs")

#x=sm2vec(file1_cor, diag = FALSE)
x=file1_cor[,1]

#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("pearson_distribution.ps", pointsize=1)

hist(x, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="freq. dist. of pears. value calcs.)", plot=TRUE)

#Other options
#breaks="Sturges" Comes up with some reasonable choice of breaks.