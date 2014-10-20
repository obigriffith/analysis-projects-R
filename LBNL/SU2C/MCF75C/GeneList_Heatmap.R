library("gplots")

#Transcript-level fpkm data from cufflinks
cuffdatafile="/Users/ogriffit/Dropbox/LBNL/Projects/SU2C/MCF75C/v2/TophatCufflinks/formattedFiles/MCF75C_fpkm.txt"

#gene-level data from Alexa-seq
alexadatafile="/Users/ogriffit/Dropbox/LBNL/Projects/SU2C/MCF75C/v2/Alexa/Matrix_GeneExpression_v53.txt"

#outfiles
outdir="/Users/ogriffit/Dropbox/LBNL/Projects/SU2C/MCF75C/v2/figures"
#heatmap_outfile = "cuff_apoptosis_heatmap.pdf"
#heatmap_outfile2 = "alexa_apoptosis_heatmap.pdf"
#heatmap_outfile = "cuff_P53signalling_heatmap.pdf"
#heatmap_outfile2 = "alexa_P53signalling_heatmap.pdf"
#heatmap_outfile = "cuff_apoptosis_heatmap2.pdf"
#heatmap_outfile2 = "alexa_apoptosis_heatmap2.pdf"
#heatmap_outfile = "cuff_E2up_E2PP2down_heatmap.pdf"
#heatmap_outfile2 = "alexa_E2up_E2PP2down_heatmap.pdf"
#heatmap_outfile = "cuff_AP1_family_heatmap.pdf"
#heatmap_outfile2 = "alexa_AP1_family_heatmap.pdf"
heatmap_outfile = "cuff_AP1_family_heatmap_noPP2.pdf"
heatmap_outfile2 = "alexa_AP1_family_heatmap_noPP2.pdf"

setwd(outdir)

#Genes of interest
#genes_interest=c("TNFRSF21","FOSL2","NUAK2","ZAK","ADAMTSL4","TP63","PMAIP1","TNFRSF11B","CXCR4","HMOX1","TGM2","FAS","LTB","TWIST1","SGK1","DNM1L","LGALS1","NTN1","BCL2L11","ATXN1","LGALS7B","HSPB8","HIPK2","CYFIP2","NGFR","PPP1R15A") #cell death
#genes_interest=c("PPP1R15A","GADD45B","ZAK","HMOX1","PMAIP1","SGK1","NUAK2","DNM1L","BCL2L11","FOSL1","FOSL2","TP63","BBC3","HIPK2","LGALS7B","LGALS1","CYFIP2","FAS","LTB","TNF","TNFRSF11B","TNFRSF21","NGFR","TWIST1","CXCR4","HSPB8","TGM2","ADAMTSL4","NTN1","ATXN1") # cell death, custom
#genes_interest=c("CDKN1A","TNFRSF10B","BBC3","BAX","SERPINE1","DDB2","TP53","MDM2","SFN","PERP","GADD45B","GADD45A") #P53 signalling
#genes_interest=c("SERPINE1","KLHDC7B","GPR87","MSN","CSF3","CYP24A1","DMP1","ACTG2","TP63","KRT17","HSPA6","FOSL1","NUAK2","AC087521.10","MAFF","TNFAIP3","FOSB","HMOX1","MGP","HSPA1B","LTB","NR4A1","CYR61","HSPA1A","LAMP3","SMTNL1","PPP1R15A","TRIM22","U4","IL28A","SNORA31","NFKBIA","IL29","IL28B","GBP1","IRF1","BATF2","ETV5","FOXI1","IFIT2","ARID5A","JUNB","HCG4P11","TNS4","IFIT3","HAS3") #E2 up, E2+PP2 down
genes_interest=c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2") #AP1 Family genes

#Experiments of interest
#experiments_interest=c("Control","E2","4OHT","E2+4OHT","PP2","E2+PP2")
experiments_interest=c("Control","E2","4OHT","E2+4OHT")

#Import data
cuffdata_import=read.table(cuffdatafile, header = TRUE, na.strings = "NA", as.is=c(1:2))
alexadata_import=read.table(alexadatafile, header = TRUE, na.strings = "NA", as.is=c(1:4))

#Break data into features info and expression values
cufffeatdata=cuffdata_import[,1:2]
cuffexpdata=cuffdata_import[,3:length(colnames(cuffdata_import))]
libnames=colnames(cuffexpdata)
libnames_clean=c("Control","E2","4OHT","E2+4OHT","PP2","E2+PP2")
colnames(cuffexpdata)=libnames_clean
rownames(cuffexpdata)=cufffeatdata[,"UCSC_Transcript_ID"]
rownames(cufffeatdata)=cufffeatdata[,"UCSC_Transcript_ID"]

alexafeatdata=alexadata_import[,1:4]
alexaexpdata=alexadata_import[,5:length(colnames(alexadata_import))]
libnames2=colnames(alexaexpdata)
libnames_clean2=c("4OHT","Control","E2+4OHT","E2+PP2","E2","PP2")
colnames(alexaexpdata)=libnames_clean2
alexaexpdata=alexaexpdata[,libnames_clean] #reorder same as for cuffdiff file
rownames(alexaexpdata)=alexafeatdata[,"EnsEMBL_Gene_ID"]
rownames(alexafeatdata)=alexafeatdata[,"EnsEMBL_Gene_ID"]

#Merge/average/add at gene symbol level
cuffexpdata_genesum=aggregate.data.frame(cuffexpdata, by=list(cufffeatdata[,"GeneSymbol"]), sum)
rownames(cuffexpdata_genesum)=cuffexpdata_genesum[,"Group.1"]
cuffexpdata_genesum=cuffexpdata_genesum[,-1]

alexaexpdata_genesum=aggregate.data.frame(alexaexpdata, by=list(alexafeatdata[,"Seq_Name"]), sum)
rownames(alexaexpdata_genesum)=alexaexpdata_genesum[,"Group.1"]
alexaexpdata_genesum=alexaexpdata_genesum[,-1]

#Extract data for genes of interest
cuffgenedata=cuffexpdata_genesum[genes_interest,]
alexagenedata=alexaexpdata_genesum[genes_interest,]

gene_names=rownames(cuffgenedata)
gene_names2=rownames(alexagenedata)

#Extract data for experiments of interest
cuffgenedata=cuffgenedata[,experiments_interest]
alexagenedata=alexagenedata[,experiments_interest]
libnames_clean=experiments_interest

#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(cuffgenedata+1)
z2 = log2(alexagenedata+1)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}


x=as.matrix(z)
main_title="Cufflinks gene-level expression"
par(cex.main=1)

pdf(file=heatmap_outfile)
#heatmap.2(t(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
#heatmap.2(t(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", Rowv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
heatmap.2(t(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
#heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=libnames_clean, labRow=gene_names, cexRow=0.9, cexCol=1, col=rev(heat.colors(75)))

dev.off()

x2=as.matrix(z2)
main_title="Alexa-seq gene-level expression"
par(cex.main=1)

pdf(file=heatmap_outfile2)
#heatmap.2(t(x2), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names2, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
#heatmap.2(t(x2), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", Rowv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names2, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
heatmap.2(t(x2), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names2, cexRow=1, cexCol=0.9, col=rev(heat.colors(75)))
#heatmap.2(x2, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=libnames_clean, labRow=gene_names2, cexRow=0.9, cexCol=1, col=rev(heat.colors(75)))

dev.off()






