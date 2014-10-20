

setwd("C:/Users/Obi/Dropbox/WashU/Projects/HCC/converge/HCC_30/CNV")
library("gplots")

cnv_genes_matrix_file="cnv.genes.matrix.tsv"

#Read in CNV data
data=read.table(cnv_genes_matrix_file, header=TRUE, sep="\t", as.is=c(1:4,7))
chromosomes=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
data=data[which(data[,"chr"]%in%chromosomes),]

data_ordered=data[mixedorder(do.call(paste, c(data[c("chr", "start")], sep = "_"))),]
matrix=data_ordered[,8:length(colnames(data))]
matrix=log2(matrix+1)

#Create fake color side bars
drug_colors=sample(c("darkorchid","darkred"), length(drug_names), replace = TRUE, prob = NULL)
subtype_colors=sample(c("red","blue","cyan","pink","yellow","green"), length(patient_ids), replace = TRUE, prob = NULL)
Mcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Ncolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Tcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
HER2colors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
PRcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
ERcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
rlab=cbind(drug_colors,drug_colors) #note, duplicated columns because function expects matrix of at least 2xn
clab=cbind(subtype_colors,Mcolors,Ncolors,Tcolors,HER2colors,PRcolors,ERcolors)
colnames(rlab)=c("Class","")
colnames(clab)=c("Subtype","M","N","T","HER2","PR","ER")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create heatmap using custom heatmap.3 source code
source("C:/Users/Obi/Dropbox/drug_predictors/Rscripts/heatmap.3.R")
pdf(file="heatmap_example.pdf")
main_title="Drug Response Predictions"
par(cex.main=1)
heatmap.3(as.matrix(matrix), na.rm = TRUE, scale="none", dendrogram="none", margins=c(4,10), Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=FALSE, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=7, KeyValueName="Value")

legend("topright",legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo"),fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()