library(hexbin)

#Read in datafiles
setwd("C:/Documents and Settings/obig/My Documents/Projects/subspace_clustering/gene_exp_distributions")

expO_Data=read.table(file="expO_exps_vs_genes.txt", header=F, sep="\t", row.names=1)
affy_Data=read.table(file="affy_exps_vs_genes.txt", header=F, sep="\t", row.names=1)
cooper_Data=read.table(file="cooper_exps_vs_genes.txt", header=F, sep="\t", row.names=1)

pdf("expO_exps_vs_genes.pdf")
hb_expO=hexbin(x=expO_Data[,2],y=expO_Data[,1])
gplot.hexbin(hb_expO, style = "nested.centroids", xlab="# genes", ylab="# experiments", lcex=0.75)
dev.off()

pdf("affy_exps_vs_genes.pdf")
hb_affy=hexbin(x=affy_Data[,2],y=affy_Data[,1])
gplot.hexbin(hb_affy, style = "nested.centroids", xlab="# genes", ylab="# experiments", lcex=0.75)
dev.off()

pdf("cooper_exps_vs_genes.pdf")
hb_cooper=hexbin(x=cooper_Data[,2],y=cooper_Data[,1])
gplot.hexbin(hb_cooper, style = "nested.centroids", xlab="# genes", ylab="# experiments", lcex=0.60)
dev.off()

#To plot three hexbin objects in different panels on the same plot
labels=vector(length=length(c(affy_Data[,1],expO_Data[,1],cooper_Data[,1])))
labels[1:length(affy_Data[,1])]="GPL96"
labels[(length(affy_Data[,1])+1):length(c(affy_Data[,1],expO_Data[,1]))]="expO"
labels[(length(c(affy_Data[,1],expO_Data[,1]))+1):length(c(affy_Data[,1],expO_Data[,1],cooper_Data[,1]))]="Cooper"

alldata=data.frame(x=c(affy_Data[,2],expO_Data[,2],cooper_Data[,2]), y=c(affy_Data[,1],expO_Data[,1],cooper_Data[,1]),a=labels)
hexbinplot(y ~ x | a, alldata, style="nested.centroids", xlab="Cluster size (# genes/cluster)", ylab="Pattern length (# dimensions/cluster)", cex.title=0.60, legend.width=0.75)



