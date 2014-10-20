#Calculate and plot correlation between Alexa-seq and Cufflinks gene-level expression estimates
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/MCF75C/v2/Alexa_vs_Cufflinks/")
outfile="Alexa_vs_Cufflinks_cors.pdf"
outfile2="Alexa_vs_Cufflinks_cors_gt1.pdf"

#Import data
cufflinks_file="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/MCF75C/v2/TophatCufflinks/formattedFiles/MCF75C_fpkm.txt"
alexaseq_file="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/MCF75C/v2/Alexa/Matrix_GeneExpression_v53.txt"
cuff_raw=read.table(file=cufflinks_file, header=TRUE, sep="\t")
alexa_raw=read.table(file=alexaseq_file, header=TRUE, sep="\t")

#Determine unique, overlapping gene set
cuff_genes_unique=as.vector(unique(cuff_raw[,"GeneSymbol"]))
alexa_genes_unique=as.vector(unique(alexa_raw[,"Seq_Name"]))
common_genes=alexa_genes_unique[which(alexa_genes_unique %in% cuff_genes_unique)]

#Extract data for just common gene set
alexa_data_common=alexa_raw[which(alexa_raw[,"Seq_Name"] %in% common_genes),c(3,4,6,9,5,7,10,8)]#Drop first two columns, and reorder
cuff_data_common=cuff_raw[which(cuff_raw[,"GeneSymbol"] %in% common_genes),]
colnames(alexa_data_common)[1]="GeneSymbol" #Give consistent column name for function to work below
lib_names=colnames(cuff_data_common)[3:8]

#Reduce redundancy, if multiple entries, take gene/transcript with highest variance, or mean across genes/transcripts?
getMeans=function(x,g,t){
 genedata=x[which(x[,"GeneSymbol"]==g),3:8]
 if(length(genedata[,1])==1){
  return(as.numeric(genedata))
 }else{
  if (t=="mean"){
   return(as.vector(mean(genedata)))
  }
  if (t=="max"){
   return(as.vector(apply(genedata,2,max)))
  }
  if (t=="sum"){
   return(as.vector(apply(genedata,2,sum)))
  }
  else{
   print("unrecognized type")
   break
  }
 }
}

test1=getMeans(cuff_data_common,"TP53","max")
test2=getMeans(alexa_data_common,"TP53","max")

alexa_data_common_mean=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
alexa_data_common_max=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
alexa_data_common_sum=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
for (i in 1:length(common_genes)){
 meandata=getMeans(alexa_data_common,common_genes[i],"mean")
 maxdata=getMeans(alexa_data_common,common_genes[i],"max")
 sumdata=getMeans(alexa_data_common,common_genes[i],"sum")
 alexa_data_common_mean[i,]=meandata
 alexa_data_common_max[i,]=maxdata
 alexa_data_common_sum[i,]=sumdata
}

cuff_data_common_mean=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
cuff_data_common_max=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
cuff_data_common_sum=matrix(NA, nrow=length(common_genes), ncol=6, dimnames=list(common_genes, lib_names))
for (i in 1:length(common_genes)){
 meandata=getMeans(cuff_data_common,common_genes[i],"mean")
 maxdata=getMeans(cuff_data_common,common_genes[i],"max")
 sumdata=getMeans(cuff_data_common,common_genes[i],"sum")
 cuff_data_common_mean[i,]=meandata
 cuff_data_common_max[i,]=maxdata
 cuff_data_common_sum[i,]=sumdata
}

#Convert values to log2, unless they are 0 in which case set them to 0
alexa_data_common_mean_log2 = log2(alexa_data_common_mean+1)
cuff_data_common_mean_log2 = log2(cuff_data_common_mean+1)
alexa_data_common_max_log2 = log2(alexa_data_common_max+1)
cuff_data_common_max_log2 = log2(cuff_data_common_max+1)
alexa_data_common_sum_log2 = log2(alexa_data_common_sum+1)
cuff_data_common_sum_log2 = log2(cuff_data_common_sum+1)


#Create plot and calculate correlation for each library for all genes
cors=vector(length=length(lib_names))
pdf(file=outfile) #all
#par (mfrow=c(3,2))
for (i in 1:length(lib_names)){
 library=lib_names[i]
 cuff_data=cuff_data_common_sum_log2[,library]
 alexa_data=alexa_data_common_sum_log2[,library]
 r=cor(x=cuff_data, y=alexa_data, method="spearman") #all
 cors[i]=r
 smoothScatter(x=cuff_data, y=alexa_data, main=library, xlab="Cufflinks", ylab="Alexa-seq")
 legend("topleft",  legend=c(paste("n =",length(common_genes)), paste("r =",format(r, digits=3))), bty="n")
}
dev.off()


#Create plot and calculate correlation for each library for all genes with value > 1 in at least one dataset
cors2=vector(length=length(lib_names))
pdf(file=outfile2) #all
#par (mfrow=c(3,2))
for (i in 1:length(lib_names)){
 library=lib_names[i]
 cuff_data=cuff_data_common_sum_log2[,library]
 alexa_data=alexa_data_common_sum_log2[,library]
 expressed_genes=which(alexa_data > 1 | cuff_data > 1)
 cuff_data_exp=cuff_data[expressed_genes]
 alexa_data_exp=alexa_data[expressed_genes]
 r=cor(x=cuff_data_exp, y=alexa_data_exp, method="spearman") #all
 cors2[i]=r
 smoothScatter(x=cuff_data_exp, y=alexa_data_exp, main=library, xlab="Cufflinks", ylab="Alexa-seq")
 legend("topleft",  legend=c(paste("n =",length(expressed_genes)), paste("r =",format(r, digits=3))), bty="n") 
}
dev.off()


