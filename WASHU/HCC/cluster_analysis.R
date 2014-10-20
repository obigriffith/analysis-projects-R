library("gplots")
library("multtest")

#Gene expression data
datafile="/Users/obigriffith/Dropbox/WashU/Projects/HCC/converge/HCC_30/RNAseq/Cufflinks_GeneLevel_FPKM.tsv"

outdir="/Users/obigriffith/Dropbox/WashU/Projects/HCC/converge/HCC_30/RNAseq/"
outfile = "HCC_gene_expression_unsupervised_heatmap.pdf"
outfile2 = "HCC_sigDE_gene_expression_log2diff_heatmap.pdf"
setwd(outdir)

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1,62:64))
header=colnames(raw_data)

exp_thresh = 1 #Minimum fpkm value to be considered expressed
pe_thresh = 0.20 #Minimum percent libraries "expressed" above exp_thresh
cov_min = 2 #Minimum coefficient of variation
cov_max = 10 #Maximum cov

#All libraries
libs=header[2:61]

#Separate individual common name from sample common name
common_names=as.data.frame(matrix(unlist(strsplit(libs,"_")), nrow=length(libs), byrow=T))

#Set up colors for sidebars
Samplecolors=as.vector(common_names[,2])
Samplecolors[Samplecolors=="normal"]="lightblue"
Samplecolors[Samplecolors=="tumor"]="blue"
clab=cbind(Samplecolors,Samplecolors)
colnames(clab)=c("","")

#Define a percent expressed function and filter out features with less than minimum
w=raw_data[,libs]
pe_fun=function(x){
 pe=length(which(x>exp_thresh))/length(x)
 return(pe)
}
pe_data=apply(w, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data_filt=raw_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
y=data_filt[,libs]
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}
cov_data=apply(y, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data_filt=data_filt[passed_cov,]

#Grab just library data for genes remaining after filtering
z=data_filt[,libs]

#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(z+1)

#Use modified heatmap.2 command to allow multiple color side bars
source("/Users/obigriffith/Dropbox/BioStars/heatmap.3.R")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmaps
#All filtered genes (unsupervised analysis)
pdf(file=outfile)
heatmap.3(as.matrix(z), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", margins=c(7,1), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, main="", labCol=libs, labRow=FALSE, cexRow=1.4, col=rev(heat.colors(75)), NumColSideColors=1, KeyValueName="log2(FPKM)")
dev.off()


#Calculate Differential expression statistics for filtered genes
#This can be done with built in stats in MTP and multiple testing correction performed at same time
sample_type=common_names[,2]
data=data_filt[,libs]
test=MTP(data,Y=sample_type,na.rm=TRUE, alternative="two.sided", test="t.pair", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP", robust=FALSE)

data_sigDE=data[which(test@adjp<0.05),]
genes_sigDE=data_filt[which(test@adjp<0.05),c(1,62:64)]
norm_data_sigDE=log2(data_sigDE[,which(sample_type=="normal")]+1)
tumor_data_sigDE=log2(data_sigDE[,which(sample_type=="tumor")]+1)
foldchanges_sigDE=tumor_data_sigDE-norm_data_sigDE

#convert_fractions=function(x){
#	x[which(x<1)]=-1/x[which(x<1)]
#	return(x)
#}
#foldchanges_converted=t(apply(foldchanges,1,convert_fractions))

#Create a custom color palatte for heatmap from yellow (down) through white (no diff) to blue (up)
nHalf=50
Min=min(foldchanges_sigDE)
Max=max(foldchanges_sigDE)
Thresh=0
rc1 <- colorRampPalette(colors = c("yellow", "white"), space="Lab")(nHalf)    
rc2 <- colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

pdf(file=outfile2)
heatmap.3(foldchanges_sigDE, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", margins=c(4,6), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, main="", labCol=as.vector(common_names[common_names[,2]=="tumor",1]), labRow=genes_sigDE[,2], cexRow=0.5, col=rampcols, breaks=rampbreaks, KeyValueName="log2 Diff")
dev.off()

data_out=cbind(genes_sigDE, test@rawp[which(test@adjp<0.05)],test@adjp[which(test@adjp<0.05)], data_sigDE)
colnames(data_out[5:6])=c("rawp","adjp")
write.table(data_out,file="HCC_FPKM_sigDE_genes.txt", sep="\t")

