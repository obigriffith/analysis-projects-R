library(geneplotter)

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/RNAseq_vs_exonarray"
setwd(outdir)
outfile="RNAseq_vs_exonarray_cors.pdf"
outfile2="RNAseq_vs_exonarray_cors_RNAseq_AB.pdf"
outfile3="RNAseq_vs_exonarray_cors_distn.pdf"
outfile4="RNAseq_vs_exonarray_cors_sample.pdf"
outfile5="RNAseq_vs_exonarray_genecors_distn.pdf"
outfile6="RNAseq_vs_exonarray_genecors_distn_50PE.pdf"
outfile7="RNAseq_vs_exonarray_genecors.txt"
outfile8="RNAseq_vs_exonarray_genecors_matrix.txt"
outfile9="RNAseq_vs_exonarray_cors_pairwise.pdf"

#RNAseq data
#gene-level
#RNAdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/RNAseq_vs_exonarray/Matrix_GeneExpression_v53.overlap.txt"
#RNAexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/RNAseq_vs_exonarray/Expressed_GeneExpression_v53.overlap.txt"

RNAdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
RNAexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"

#exonarray data
#gene-level
#EAdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/RNAseq_vs_exonarray/breastExon_genelevel.overlap.csv"
EAdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ExonArray/prefiltered/breastExon_genelevel.csv"

RNA_data_import=read.table(RNAdatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4), check.names=FALSE)
RNA_exp_status_import=read.table(RNAexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4), check.names=FALSE)
EA_data_import=read.csv(EAdatafile, check.names=FALSE) 

#Break RNAseq data into features info and expression values
RNA_feat_data=RNA_data_import[,1:4]
RNA_data=RNA_data_import[,5:length(colnames(RNA_data_import))]
RNA_exp_status=RNA_exp_status_import[,5:length(colnames(RNA_data_import))]

#Fix misnamed libraries in RNA_data and RNA_exp_status
colnames(RNA_data)[which(colnames(RNA_data)=="MDAMB13v1")]="MDAMB134VI"
colnames(RNA_data)[which(colnames(RNA_data)=="SUM1315")]="SUM1315MO2"
colnames(RNA_exp_status)=colnames(RNA_data)
rownames(RNA_data)=RNA_feat_data[,4] #Use ENSG as rowname because this is unique
rownames(RNA_exp_status)=RNA_feat_data[,4] #Use ENSG as rowname because this is unique
rownames(RNA_feat_data)=RNA_feat_data[,4] #Use ENSG as rowname because this is unique

#separate gene names and data for exon array data
EA_feat_data=EA_data_import[,1]
EA_data=EA_data_import[,2:length(colnames(EA_data_import))]
rownames(EA_data)=EA_feat_data

#Determine overlapping libraries/cell lines between RNA and exonArray data
RNAseq_exonarray_libs=sort(colnames(EA_data)[colnames(EA_data)%in%colnames(RNA_data)])

#Determine overlapping genes
#First need to simplify RNAseq data to eliminate redundant gene symbols (choose most variant)
#Calculate COV for all genes
cov_fun=function(x){
  cov=sd(x, na.rm=TRUE)/abs(mean(x, na.rm=TRUE))
  return(cov)
}
RNA_data_cov=apply(RNA_data,1,cov_fun)

#Combine with RNA feature data
RNA_data_cov=cbind(RNA_feat_data,RNA_data_cov)

#Get unique gene symbols
RNA_uniq_symbols=unique(RNA_feat_data[,"Seq_Name"])

#For each unique gene symbol get the gene with highest COV
max_cov_gene_fun=function(x){
 gene_data=RNA_data_cov[RNA_data_cov[,"Seq_Name"]==x,]
 max_cov_gene=gene_data[order(gene_data[,"RNA_data_cov"], decreasing=TRUE)[1],"EnsEMBL_Gene_ID"]
 return(max_cov_gene)
}
max_cov_genes=sapply(RNA_uniq_symbols, max_cov_gene_fun)

#Filter down RNA data to max cov gene set
RNA_feat_data=RNA_feat_data[max_cov_genes,]
RNA_data=RNA_data[max_cov_genes,]
RNA_exp_status=RNA_exp_status[max_cov_genes,]

#Reset rownames to gene symbol which should now be unique
rownames(RNA_data)=RNA_feat_data[,3] #Use ENSG as rowname because this is unique
rownames(RNA_exp_status)=RNA_feat_data[,3] #Use ENSG as rowname because this is unique
rownames(RNA_feat_data)=RNA_feat_data[,3] #Use ENSG as rowname because this is unique

#determine overlap between RNA max cov genes and exon array data
RNA_EA_overlap_genes=as.vector(EA_feat_data[EA_feat_data%in%RNA_feat_data[,3]])

#Filter down data to just overlap genes
EA_data=EA_data[RNA_EA_overlap_genes,]
RNA_data=RNA_data[RNA_EA_overlap_genes,]
RNA_feat_data=RNA_feat_data[RNA_EA_overlap_genes,]
RNA_exp_status=RNA_exp_status[RNA_EA_overlap_genes,]

#Create vector of gene names
genes=as.vector(RNA_feat_data[,3])

#RNAseq data needs to be log2'd?
RNA_data=log2(RNA_data+1)

#Create plot and calculate correlation for each library for all genes
cors=vector(length=length(RNAseq_exonarray_libs))
pdf(file=outfile) #all
par (mfrow=c(3,3))
for (i in 1:length(RNAseq_exonarray_libs)){
 library=RNAseq_exonarray_libs[i]
 exonArray_data=EA_data[,library]
 RNAseq_data=RNA_data[,library]
 r=cor(x=RNAseq_data, y=exonArray_data, method="spearman") #all
 cors[i]=r
 smoothScatter(x=RNAseq_data, y=exonArray_data, main=library, xlab="RNA-seq", ylab="Exon Array")
 legend("topleft",  legend=paste("r=",format(r, digits=3)), bty="n") 
}
dev.off()

#####
#NOTE: Would be good to update to also filter out genes not expressed above background in exon-array data
#Problem, how to define background in exon array data?
#####


#Create plot and calculate correlation for each library for all genes with RNAseq expression above background
cors2=vector(length=length(RNAseq_exonarray_libs))
pdf(file=outfile2) #RNAseq expressed
par (mfrow=c(3,3))
for (i in 1:length(RNAseq_exonarray_libs)){
 library=RNAseq_exonarray_libs[i]
 exonArray_data=EA_data[,library]
 RNAseq_data=RNA_data[,library]
 RNAseq_expr_data=RNA_exp_status[,library]
 RNAseq_expressed=which(RNAseq_expr_data==1)
 r=cor(x=RNAseq_data[RNAseq_expressed], y=exonArray_data[RNAseq_expressed], method="spearman") #RNAseq expressed
 cors2[i]=r
 smoothScatter(x=RNAseq_data[RNAseq_expressed], y=exonArray_data[RNAseq_expressed], main=library, xlab="RNA-seq", ylab="Exon Array") #RNAseq expressed
 legend("topleft",  legend=c(paste("r=",format(r, digits=3), sep=""), paste("n=",length(RNAseq_data[RNAseq_expressed]), sep="")), bty="n") #RNAseq expressed
}
dev.off()

#Calculate correlation between each library and all other libraries to make sure that library matches itself best
cors3=matrix(data=NA, nrow=length(RNAseq_exonarray_libs), ncol=length(RNAseq_exonarray_libs), byrow=FALSE, dimnames=NULL)
rownames(cors3)=RNAseq_exonarray_libs
colnames(cors3)=RNAseq_exonarray_libs
pdf(file=outfile9)
par (mfrow=c(2,2))
for (i in 1:length(RNAseq_exonarray_libs)){
 libraryA=RNAseq_exonarray_libs[i]
 RNAseq_data=RNA_data[,libraryA]
 RNAseq_expr_data=RNA_exp_status[,libraryA]
 RNAseq_expressed=which(RNAseq_expr_data==1)
 for (j in 1:length(RNAseq_exonarray_libs)){
  libraryB=RNAseq_exonarray_libs[j]
  exonArray_data=EA_data[,libraryB]
  r=cor(x=RNAseq_data[RNAseq_expressed], y=exonArray_data[RNAseq_expressed], method="spearman") #all
  cors3[i,j]=r
 }
 lib_cors_sorted=sort(cors3[i,], decreasing=TRUE)
 plot(lib_cors_sorted, xaxt = "n", xlab=NA, ylab="Spearman rho", main=libraryA)
 axis(1, at=c(1:length(lib_cors_sorted)), labels=names(lib_cors_sorted), las=3, cex.axis=0.5, xlab=FALSE) 
}
dev.off()
write.table(cors3, file=outfile8, sep="\t")

#Create histogram of distn of rho values
pdf(file=outfile3)
hist(cors2, xlim=c(-1,1), col="blue", main="Distribution of correlations for RNAseq vs exon array", xlab="Spearman rho")
legend_text=c(
paste("mean rho =",format(mean(cors2, na.rm=TRUE), digits=3), sep=" "),
paste("range:",format(min(cors2), digits=3),"to",format(max(cors2), digits=3), sep=" ")
)
legend("topleft", legend=legend_text, bty="n")
dev.off()

#Create individual figure for example with mean rho value
pdf(file=outfile4)
library="HCC1500"
exonArray_data=EA_data[,library]
RNAseq_data=RNA_data[,library]
RNAseq_expr_data=RNA_exp_status[,library]
RNAseq_expressed=which(RNAseq_expr_data==1)
r=cor(x=RNAseq_data, y=exonArray_data, method="spearman") #all
cors[i]=r
smoothScatter(x=RNAseq_data, y=exonArray_data, main=library, xlab="RNA-seq", ylab="Exon Array")
legend("topleft",  legend=paste("r=",format(r, digits=3)), bty="n") 
dev.off()


#Calculate correlations on a gene-by-gene basis for all libraries (with both datatypes)
#First, trim down to just libraries shared between datatypes
RNA_data_overlap=RNA_data[,RNAseq_exonarray_libs]
EA_data_overlap=EA_data[,RNAseq_exonarray_libs]
RNA_exp_status_overlap=RNA_exp_status[,RNAseq_exonarray_libs]

#Create a function to calculate correlation for a single gene across all libraries
genecor_fun=function(x){ 
 r=cor(x=as.numeric(RNA_data_overlap[x,]), y=as.numeric(EA_data_overlap[x,]), method="spearman")
 return(r)
}
gene_cors=apply(as.array(genes), 1, genecor_fun)

pdf(file=outfile5)
hist(gene_cors, col="blue", xlim=c(-1,1), main="Distribution of gene correlations across all cell lines", xlab="Spearman rho")
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
pe_data=apply(RNA_exp_status_overlap, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
RNA_data_overlap_filt=RNA_data_overlap[passed_pe,]
EA_data_overlap_filt=EA_data_overlap[passed_pe,]
genes_filt=rownames(RNA_data_overlap_filt)
gene_cors_filt=apply(as.array(genes_filt), 1, genecor_fun)

pdf(file=outfile6)
hist(gene_cors_filt, col="blue", xlim=c(-1,1), main="Distribution of gene correlations across all cell lines", xlab="Spearman rho")
legend_text=c(
paste("genes expressed in >",pe_thresh*100,"% cell lines (RNAseq)", sep=""),
paste("n =",length(genes_filt),"genes total", sep=" "),
paste("mean rho =",format(mean(gene_cors_filt), digits=3), sep=" "),
paste(format(length(which(gene_cors_filt>0.5))/length(gene_cors_filt)*100, digits=3), "% genes with rho > 0.5", sep="")
)
legend("topleft", legend=legend_text, bty="n")
dev.off()

#Write results to file:
data=cbind(genes,gene_cors,pe_data)
colnames(data)=c("Gene","rho","% cell lines expressed (RNAseq)")
write.table(data, file=outfile7, row.names=FALSE, sep="\t")

