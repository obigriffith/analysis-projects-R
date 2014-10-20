library("gplots")
library("heatmap.plus")

#Gene expression data
datafile="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/RNAseq_heatmap_analysis/Cufflinks_GeneLevel_Fixed.tsv"

outdir="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/RNAseq_heatmap_analysis/"
outfile = "BRAF_gene_expression_unsupervised_heatmap.pdf"
setwd(outdir)

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1,12:14))
header=colnames(raw_data)

exp_thresh = 10 #Minimum fpkm value to be considered expressed
pe_thresh = 0.10 #Minimum percent libraries "expressed" - at least 2/10 libraries
cov_min = 1 #Minimum coefficient of variation
cov_max = 10 #Maximum cov

#All libraries
libs=header[2:11]
lib_ids=libs
lib_names=c("F6-8","F6-8/R12","F6-8/R6","F6-8/R8","Vaco432","Vaco432/R10","Vaco432/R19","Vaco432/R22","Vaco432/R23","Vaco432/R1")
lib_types=c("F6Parent","F6","F6","F6","VacoParent","Vaco","Vaco","Vaco","Vaco","Vaco")

#Set up colors for sidebars
Samplecolors=lib_types
Samplecolors[Samplecolors=="VacoParent"]="lightblue"
Samplecolors[Samplecolors=="Vaco"]="blue"
Samplecolors[Samplecolors=="F6Parent"]="lightgreen"
Samplecolors[Samplecolors=="F6"]="green"
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
source("/Users/ogriffit/Dropbox/drug_predictors/Rscripts/heatmap.3.R")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Manually set order for this heatmap
col_order=c(4,1,2,3,10,7,6,9,5,8)
coldd=as.dendrogram(myclust(mydist(as.matrix(t(z)))))
coldd=reorder(coldd,col_order)

#Heatmaps
#All filtered genes (unsupervised analysis)
pdf(file=outfile)
heatmap.3(as.matrix(z), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", margins=c(7,1), Rowv=TRUE, Colv=coldd, ColSideColors=clab, symbreaks=FALSE, main="", labCol=lib_names, labRow=FALSE, cexRow=1.4, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="log2(FPKM)")
dev.off()




