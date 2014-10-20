library(geneplotter)

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/Pancreas/DSN_RNAseq_vs_U133plus2/customCDF"
setwd(outdir)

outfile="RNAseq_vs_U133plus2_cors.pdf"
outfile2="RNAseq_vs_U133plus2_cors_RNAseq_AB.pdf"
outfile3="RNAseq_vs_U133plus2_cors_distn.pdf"
outfile4="RNAseq_vs_U133plus2_cors_sample.pdf"
outfile5="RNAseq_vs_U133plus2_genecors_distn.pdf"
outfile6="RNAseq_vs_U133plus2_genecors_distn_50PE.pdf"
outfile7="RNAseq_vs_U133plus2_genecors.txt"

RNAseq_libs=c("HCT20639","HCT20640","HCT20641","HCT20642")
U133APlus2_libs=c("AGR06.1126.CEL","AGR06.1127.CEL","AGR06.1128.CEL","AGR06.1129.CEL")

#RNAseq data
#gene-level
RNAdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/Pancreas/DSN_RNAseq_vs_U133plus2/Matrix_GeneExpression_v53.txt"
RNAexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/Pancreas/DSN_RNAseq_vs_U133plus2/Expressed_GeneExpression_v53.txt"

#exonarray data
#gene-level
U133Plus2datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Pancreas/DSN_RNAseq_vs_U133plus2/ALL_gcrma_customCDF.txt"

#Import data
RNA_data_import=read.table(RNAdatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
RNA_exp_status_import=read.table(RNAexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
U133Plus2_data_import=read.table(U133Plus2datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))

#Break RNAseq data into features info and expression values
RNA_feat_data=RNA_data_import[,1:4]
RNA_data=RNA_data_import[,5:length(colnames(RNA_data_import))]
RNA_exp_status=RNA_exp_status_import[,5:length(colnames(RNA_data_import))]
rownames(RNA_data)=RNA_feat_data[,4]
rownames(RNA_exp_status)=RNA_feat_data[,4]
rownames(RNA_feat_data)=RNA_feat_data[,4]


#Break array data into features info and expression values
U133Plus2_feat_data=U133Plus2_data_import[,1:3]
U133Plus2_data=U133Plus2_data_import[,4:length(colnames(U133Plus2_data_import))]
rownames(U133Plus2_data)=U133Plus2_feat_data[,1]
rownames(U133Plus2_feat_data)=U133Plus2_feat_data[,1]

#find unique, overlapping set of genes
overlap_genes=unique(U133Plus2_feat_data[U133Plus2_feat_data[,3] %in% RNA_feat_data[,3],3])

#Need to choose just one array probe for each gene (use most variable)
max_var_probes=rep(NA,length(overlap_genes))
for (i in 1:length(overlap_genes)){
 gene=overlap_genes[i]
 probes=U133Plus2_feat_data[which(U133Plus2_feat_data[,3]==gene),1]
 probe_data=U133Plus2_data[probes,]
 probe_vars=apply(probe_data, 1, var)
 max_var_probe=probes[order(probe_vars, decreasing=TRUE)[1]]
 max_var_probes[i]=max_var_probe
}

#Filter array data down to just max_var probes and samples of interest
U133Plus2_data=U133Plus2_data[max_var_probes,U133APlus2_libs]
U133Plus2_feat_data=U133Plus2_feat_data[max_var_probes,]

#Set library (colnames) to same as for RNAseq data
colnames(U133Plus2_data)=RNAseq_libs

#Need to choose just one ENSG for each gene symbol (use most variable)
max_var_ENSGs=rep(NA,length(overlap_genes))
for (i in 1:length(overlap_genes)){
 gene=overlap_genes[i]
 ENSGs=RNA_feat_data[which(RNA_feat_data[,3]==gene),4]
 ENSG_data=RNA_data[ENSGs,]
 ENSG_vars=apply(ENSG_data, 1, var)
 max_var_ENSG=ENSGs[order(ENSG_vars, decreasing=TRUE)[1]]
 max_var_ENSGs[i]=max_var_ENSG
}

#Filter RNAseq data down to just max_var ENSGs
RNA_feat_data=RNA_feat_data[max_var_ENSGs,]
RNA_data=RNA_data[max_var_ENSGs,]
RNA_exp_status=RNA_exp_status[max_var_ENSGs,]

#Create vector of gene names
genes=as.vector(RNA_feat_data[,3])

#Now that both datasets are filtered down to same set of unique genes, use gene symbols as rownames
rownames(RNA_data)=genes
rownames(U133Plus2_data)=genes

#RNAseq data needs to be log2'd?
RNA_data=log2(RNA_data+1)

#Create plot and calculate correlation for each library for all genes
cors=vector(length=length(RNAseq_libs))
pdf(file=outfile) #all
par (mfrow=c(2,2))
for (i in 1:length(RNAseq_libs)){
 library=RNAseq_libs[i]
 Array_data=U133Plus2_data[,library]
 RNAseq_data=RNA_data[,library]
 r=cor(x=RNAseq_data, y=Array_data, method="spearman") #all
 cors[i]=r
 smoothScatter(x=RNAseq_data, y=Array_data, main=library, xlab="RNA-seq", ylab="U133Plus2")
 legend("topleft",  legend=paste("r=",format(r, digits=3)), bty="n") 
}
dev.off()

#Create plot and calculate correlation for each library for all genes with RNAseq expression above background
cors2=vector(length=length(RNAseq_libs))
pdf(file=outfile2) #RNAseq expressed
par (mfrow=c(2,2))
for (i in 1:length(RNAseq_libs)){
 library=RNAseq_libs[i]
 Array_data=U133Plus2_data[,library]
 RNAseq_data=RNA_data[,library]
 RNAseq_expr_data=RNA_exp_status[,library]
 RNAseq_expressed=which(RNAseq_expr_data==1)
 r=cor(x=RNAseq_data[RNAseq_expressed], y=Array_data[RNAseq_expressed], method="spearman") #RNAseq expressed
 cors2[i]=r
 smoothScatter(x=RNAseq_data[RNAseq_expressed], y=Array_data[RNAseq_expressed], main=library, xlab="RNA-seq", ylab="U133Plus2") #RNAseq expressed
 legend("topleft",  legend=c(paste("r=",format(r, digits=3), sep=""), paste("n=",length(RNAseq_data[RNAseq_expressed])," expressed (RNA-seq)", sep="")), bty="n") #RNAseq expressed
}
dev.off()

#Calculate correlations on a gene-by-gene basis for all libraries (with both datatypes)
#Create a function to calculate correlation for a single gene across all libraries
genecor_fun=function(i){ 
 r=cor(x=as.numeric(RNA_data[i,]), y=as.numeric(U133Plus2_data[i,]), method="spearman")
 return(r)
}
gene_cors=apply(as.array(genes), 1, genecor_fun)


pdf(file=outfile5)
hist(gene_cors, col="blue", xlim=c(-1,1), main="Distribution of gene correlations across 4 tumors", xlab="Spearman rho")
legend_text=c(
paste("n =",length(genes),"genes total", sep=" "),
paste("mean rho =",format(mean(gene_cors, na.rm=TRUE), digits=3), sep=" "),
paste(format(length(which(gene_cors>0.5))/length(gene_cors)*100, digits=3), "% genes with rho > 0.5", sep="")
)
legend("topleft", legend=legend_text, bty="n")
dev.off()



#Look at only genes with at least X% expressed above background (from RNAseq)
#Define a percent expressed function
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}

#apply filter on percent expressed (RNAseq data)
pe_thresh=0.5
pe_data=apply(RNA_exp_status, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
RNA_data_filt=RNA_data[passed_pe,]
U133Plus2_data_filt=U133Plus2_data[passed_pe,]
genes_filt=rownames(RNA_data_filt)
gene_cors_filt=apply(as.array(genes_filt), 1, genecor_fun)


pdf(file=outfile6)
hist(gene_cors_filt, col="blue", xlim=c(-1,1), main="Distribution of gene correlations across all 4 tumors", xlab="Spearman rho")
legend_text=c(
paste("genes expressed in >",pe_thresh*100,"% cell lines (RNAseq)", sep=""),
paste("n =",length(genes_filt),"genes total", sep=" "),
paste("mean rho =",format(mean(gene_cors_filt, na.rm=TRUE), digits=3), sep=" "),
paste(format(length(which(gene_cors_filt>0.5))/length(gene_cors_filt)*100, digits=3), "% genes with rho > 0.5", sep="")
)
legend("topleft", legend=legend_text, bty="n")
dev.off()

#Write results to file:
data=cbind(genes,gene_cors,pe_data)
colnames(data)=c("Gene","rho","% cell lines expressed (RNAseq)")
write.table(data, file=outfile7, row.names=FALSE, sep="\t")






#Unused code
#Create individual figure for example with mean rho value
pdf(file=outfile4)
library="HCC1500"
exonArray_data=U133Plus2_data[,library]
RNAseq_data=RNA_data[,library]
RNAseq_expr_data=RNA_exp_status[,library]
RNAseq_expressed=which(RNAseq_expr_data==1)
r=cor(x=RNAseq_data, y=exonArray_data, method="spearman") #all
cors[i]=r
smoothScatter(x=RNAseq_data, y=exonArray_data, main=library, xlab="RNA-seq", ylab="Exon Array")
legend("topleft",  legend=paste("r=",format(r, digits=3)), bty="n") 
dev.off()


