library("gplots")

#gene-level data from Alexa-seq
alexadatafile="/Users/ogriffit/Dropbox/LBNL/Projects/SU2C/MCF7/Alexa/Matrix_GeneExpression_v53.txt"

#outfiles
outdir="/Users/ogriffit/Dropbox/LBNL/Projects/SU2C/MCF7/figures"
heatmap_outfile = "MCF7_WS8_alexa_heatmap.pdf"

setwd(outdir)

#Genes of interest
#genes_interest=c("RHOF", "EHD2", "KIF12", "FKBP4", "SLC25A25", "SYTL5", "CHST8", "STC2", "TGM2", "FRK", "FHL2", "HOMER3", "PHLDA2", "RDX", "RASIP1", "CDC42EP2", "TUBB6", "PTAFR", "DCBLD2", "SLC16A1", "KCTD12", "GAL", "LRP8", "ASPH", "SLMAP", "MESP1", "C3", "MRPS10", "GTPBP3", "ALDH1L2", "PTP4A1", "SH2D5", "LEPREL1", "AIDA", "ENO2", "ACTL8", "LOXL3", "ITGB1")
genes_interest=c("RHOF", "EHD2", "KIF12", "SLC25A25", "SYTL5", "CHST8", "STC2", "TGM2", "FRK", "FHL2", "HOMER3", "PHLDA2", "RDX", "RASIP1", "CDC42EP2", "TUBB6", "PTAFR", "DCBLD2", "SLC16A1", "KCTD12", "GAL", "LRP8", "ASPH", "SLMAP", "MESP1", "C3", "MRPS10", "GTPBP3", "ALDH1L2", "PTP4A1", "SH2D5", "LEPREL1", "AIDA", "ENO2", "ACTL8", "LOXL3", "ITGB1")


#Experiments of interest
experiments_interest=c("WS8_Control","WS8_E2","WS8_4OHT") #WS8 data

#Import data
alexadata_import=read.table(alexadatafile, header = TRUE, na.strings = "NA", as.is=c(1:4))

#Break data into features info and expression values
alexafeatdata=alexadata_import[,1:4]
alexaexpdata=alexadata_import[,5:length(colnames(alexadata_import))]
libnames=colnames(alexaexpdata)

#CHECK ORDER OF CLEAN NAMES AGAINST IMPORTED DATA ABOVE
libnames_clean=c("5C_4OHT","5C_Control","5CE2_4OHT","5C_E2","5CE2_Control","5CE2_E2","5CE2PP2_4OHT","5CE2PP2_Control","5CE2PP2_E2","5CPP2_4OHT","5CPP2_Control","5CPP2_E2","WS8_4OHT","WS8_Control","WS8_E2")

colnames(alexaexpdata)=libnames_clean
rownames(alexaexpdata)=alexafeatdata[,"EnsEMBL_Gene_ID"]
rownames(alexafeatdata)=alexafeatdata[,"EnsEMBL_Gene_ID"]

#Merge/average/add at gene symbol level
alexaexpdata_genesum=aggregate.data.frame(alexaexpdata, by=list(alexafeatdata[,"Seq_Name"]), sum)
rownames(alexaexpdata_genesum)=alexaexpdata_genesum[,"Group.1"]
alexaexpdata_genesum=alexaexpdata_genesum[,-1]

#Extract data for genes of interest
alexagenedata=alexaexpdata_genesum[genes_interest,]
gene_names=rownames(alexagenedata)

#Extract data for experiments of interest
alexagenedata=alexagenedata[,experiments_interest]
libnames_clean=experiments_interest

#Convert values to log2, unless they are 0 in which case set them to 0
z=log2(alexagenedata+1)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

x=as.matrix(z)
#main_title="Alexa-seq gene-level expression"
main_title=""
par(cex.main=1)

experiment_display_names=c("Control","E2","4-OHT")
pdf(file=heatmap_outfile)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm=TRUE, scale="none", dendrogram="none", Rowv=TRUE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=gene_names, labCol=experiment_display_names, cexRow=0.9, cexCol=1.3, col=rev(heat.colors(75)), margins = c(5, 20), keysize=1.6)
dev.off()






