library(mclust)
library("heatmap.plus")
library(gplots)

datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/processing/processed_final2/test_train_survival/combined/ALL_gcrma.txt"
#datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/processing/processed_final2/test_train_survival/customCDF/ALL_gcrma.txt"
clindatafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.txt"
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/single_gene/ZNF217"
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/single_gene/ZNF217/update/"

outfile="test_train_survival_ALL_gcrma_ZNF217_w_clin_data.txt"
outfile2="test_train_survival_ALL_gcrma_ZNF217_mix_model.pdf"
outfile3="test_train_survival_ALL_gcrma_ZNF217_hist.pdf"
outfile4="test_train_survival_ALL_gcrma_genesInterest_heatmap.pdf"
outfile5="test_train_survival_ALL_gcrma_genesInterest_cors.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clin_data=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
raw_data=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

#Change to output dir
setwd(outdir)
header=colnames(raw_data)

#Define function to get single-gene Data
getGeneData=function(x) {
 Gene_Data_all=raw_data[which(raw_data[,3]==x),]
 Gene_exp_data=Gene_Data_all[,4:length(header)]
 probe_names=Gene_Data_all[,1]
 gene_names=Gene_Data_all[,3]
 names=paste(probe_names,"(",gene_names,")",sep="")
 rownames(Gene_exp_data)=names
 return(Gene_exp_data)
}

#Get gene data for ZNF217
ZNF217_exp_data=getGeneData("ZNF217")

#Other genes of interest to correlate with ZNF217
#E-cadherin, Twist, Snail1, Snail2, Twist2?, Vimentin, Keratin 8, Keratin 18, Keratin 17, Keratin 19, Keratin 5, Keratin 7, Keratin 14, ESR1
#CDH1, TWIST1, SNAI1, SNAI2, TWIST2, VIM, KRT8, KRT18, KRT17, KRT19, KRT5, KRT7, KRT14, ESR1

CDH1_exp_data=getGeneData("CDH1")
TWIST1_exp_data=getGeneData("TWIST1")
SNAI1_exp_data=getGeneData("SNAI1")
SNAI2_exp_data=getGeneData("SNAI2")
VIM_exp_data=getGeneData("VIM")
KRT8_exp_data=getGeneData("KRT8")
KRT18_exp_data=getGeneData("KRT18")
KRT17_exp_data=getGeneData("KRT17")
KRT19_exp_data=getGeneData("KRT19")
KRT5_exp_data=getGeneData("KRT5")
KRT7_exp_data=getGeneData("KRT7")
KRT14_exp_data=getGeneData("KRT14")
ESR1_exp_data=getGeneData("ESR1")
ERBB3_exp_data=getGeneData("ERBB3")

Gene_exp_data=rbind(ZNF217_exp_data,CDH1_exp_data,TWIST1_exp_data,SNAI1_exp_data,SNAI2_exp_data,VIM_exp_data,KRT8_exp_data,KRT18_exp_data,KRT17_exp_data,KRT19_exp_data,KRT5_exp_data,KRT7_exp_data,KRT14_exp_data,ESR1_exp_data,ERBB3_exp_data)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap (single color sidebar) - ERBB2/GRB7/ESR1 only
x=as.matrix(Gene_exp_data)
row_names=rownames(x)

pdf(file=outfile4)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="column", Rowv=FALSE, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#From heatmap chose best/representative probe for each gene
Gene_exp_data_best=Gene_exp_data[c("7764_at(ZNF217)","201131_s_at(CDH1)","7291_at(TWIST1)","6615_at(SNAI1)","6591_at(SNAI2)","7431_at(VIM)","209008_x_at(KRT8)","3875_at(KRT18)","3872_at(KRT17)","3880_at(KRT19)","3852_at(KRT5)","3855_at(KRT7)","3861_at(KRT14)","205225_at(ESR1)","2065_at(ERBB3)"),]

#Calculate correlations between ZNF217 and other genes of interest
#Get list of genes and loop through
probes=rownames(Gene_exp_data_best)

#Create dataframe to store results
cors = data.frame(cbind(probes,r_spearman=NA,p_spearman=NA,r_pearson=NA,p_pearson=NA), stringsAsFactors=FALSE)

for (i in 1:length(probes)){
  probe=probes[i]
  cors[i,"r_spearman"]=cor(x=as.numeric(Gene_exp_data_best["7764_at(ZNF217)",]), y=as.numeric(Gene_exp_data_best[probe,]), method="spearman")
  cors[i,"p_spearman"]=cor.test(x=as.numeric(Gene_exp_data_best["7764_at(ZNF217)",]), y=as.numeric(Gene_exp_data_best[probe,]), method="spearman", alternative = "two.sided")$p.value
  cors[i,"r_pearson"]=cor(x=as.numeric(Gene_exp_data_best["7764_at(ZNF217)",]), y=as.numeric(Gene_exp_data_best[probe,]), method="pearson")
  cors[i,"p_pearson"]=cor.test(x=as.numeric(Gene_exp_data_best["7764_at(ZNF217)",]), y=as.numeric(Gene_exp_data_best[probe,]), method="pearson", alternative = "two.sided")$p.value
}

write.table(cors, file=outfile5, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#Attempt to separate ZNF217 data based on gaussian mixture model
#Allow mclust to determine best model - Answer in this case was 2 clusters (one in middle and one with both extremes). Break into 3.
x=as.numeric(ZNF217_exp_data[2,])
mclust_ZNF217=Mclust(x)
summary(mclust_ZNF217, x) #Gives you list of values returned
classification_ZNF217=mclust_ZNF217$classification
num_clusters=mclust_ZNF217$G

pdf(file=outfile2)
par(mfrow=c(2,1), oma=c(2,2,2,2))
hist(x[classification_ZNF217==1], xlim=c(2,15), col="blue", xlab="Log2 GCRMA value", main="ZNF217 components to retain")
hist(x[classification_ZNF217==2], xlim=c(2,15), col="red", xlab="Log2 GCRMA value", main="ZNF217 components to filter")
title(main="Separation of ZNF217 by model-based clustering", outer=TRUE)
dev.off()

#Choose cutoffs - use max and min of "middle" cluster (#2) to break into three components
ZNF217_cutoff1=min(x[classification_ZNF217==2])
ZNF217_cutoff2=max(x[classification_ZNF217==2])
ZNF217_low_count=length(which(x <= ZNF217_cutoff1))
ZNF217_int_count=length(which(x > ZNF217_cutoff1 & x < ZNF217_cutoff2))
ZNF217_high_count=length(which(x >= ZNF217_cutoff2))

#Histogram - ZNF217
pdf(file=outfile3)
hist(as.numeric(x), xlim=c(2,15), breaks=40, col="blue", main=paste("ZNF217 (7764_at) histogram for all 858 samples", sep=""), xlab="Log2 GCRMA value", ylab="Frequency")
abline(v=ZNF217_cutoff1, col="red")
abline(v=ZNF217_cutoff2, col="red")
legend("left", legend=c(paste("low cutoff=",ZNF217_cutoff1,sep=""),paste("high cutoff=",ZNF217_cutoff2,sep=""),paste("N_low=",ZNF217_low_count,sep=""),paste("N_int=",ZNF217_int_count,sep=""),paste("N_high=",ZNF217_high_count,sep="")),bty="n")
dev.off()

#Extract low, int, high samples
ZNF217_low=which(x <= ZNF217_cutoff1)
ZNF217_int=which(x > ZNF217_cutoff1 & x < ZNF217_cutoff2)
ZNF217_high=which(x >= ZNF217_cutoff2)

#Join clinical and gene data and write to file
clin_gene_data=cbind(clin_data,t(ZNF217_exp_data))

#Create new column to store cutoff groups 
clin_gene_data[ZNF217_low,"ZNF217_group"]="low"
clin_gene_data[ZNF217_int,"ZNF217_group"]="int"
clin_gene_data[ZNF217_high,"ZNF217_group"]="high"

write.table(clin_gene_data, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




