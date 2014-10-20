library(hexbin)

#Read in datafiles
setwd("C:/Documents and Settings/obig/My Documents/Projects/subspace_clustering/StanfordPromoters/NegControlAnalysis/NegOnly_vs_PosOnly")
NegData=read.table(file="cooper.expander.shortNames.probeId.20060328.NEG_genes_vs_dims.txt", header=T, sep="\t", row.names=1)
PosData=read.table(file="cooper.expander.shortNames.probeId.20060328.POS_genes_vs_dims.txt", header=T, sep="\t", row.names=1)

#To Summarize numbers at each gene/dim level
NegDataTable=table(NegData)
PosDataTable=table(PosData)

#Get x and y data for two datasets
#Negative data
neg_genes=NegData[,1]
neg_dims=NegData[,2]
#Positive data
pos_genes=PosData[,1]
pos_dims=PosData[,2]

#Create hexbin objects for two data distributions
hb_neg=hexbin(x=neg_genes, y=neg_dims)
hb_pos=hexbin(x=pos_genes, y=pos_dims)

#Plot two hexbin objects on different plots
plot(hb_neg, colramp=function(n){BTY(n,beg=15,end=225)})
plot(hb_pos, colramp=function(n){magent(n,beg=15,end=225)})

#To plot two hexbin objects in different panels on the same plot
labels=vector(length=length(c(neg_genes,pos_genes)))
labels[1:length(neg_genes)]="Negatives"
labels[length(neg_genes)+1:length(pos_genes)]="Positives"

pdf(file="Negs_vs_Positives_hexbin.pdf")
alldata=data.frame(x=c(neg_genes,pos_genes), y=c(neg_dims,pos_dims),a=labels)
hexbinplot(y ~ x | a, alldata, style="nested.centroids", xlab="Cluster size (# genes/cluster)", ylab="Pattern length (# dimensions/cluster)", cex.title=0.60, legend.width=0.75)
dev.off()

#To plot the difference between two hexbin objects on the same plot
#Set the x and y boundaries for binning and plotting
xbnds <- range(c(neg_genes,pos_genes), na.rm = TRUE)
ybnds <- range(c(neg_dims,pos_dims), na.rm = TRUE)
plot.new()
hb_neg=hexbin(x=neg_genes, y=neg_dims, xbnds=xbnds, ybnds=ybnds)
hb_pos=hexbin(x=pos_genes, y=pos_dims, xbnds=xbnds, ybnds=ybnds)
hdiffplot(erode(hb_neg, cdfcut = 0.25), erode(hb_pos, cdfcut = 0.25),unzoom = 1.3)

