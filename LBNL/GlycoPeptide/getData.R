outdir="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide"
setwd(outdir)

#genes of interest
gene_interest_file="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide/HumanGlycosyltransferases_etc.txt"
gene_interest_import=read.table(gene_interest_file, header = TRUE, na.strings = "NA", sep="\t", as.is=1:2)
genes_interest=gene_interest_import[,1]
rownames(gene_interest_import)=gene_interest_import[,1]

#Cell line info
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/Projects/GlycoPeptide/BCCL_data_list.txt"
cell_data_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=1)
lines=cell_data_import[,1]

### exonarray data ###
EAdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ExonArray/prefiltered/breastExon_genelevel.csv"
EA_data_import=read.csv(EAdatafile) 
EA_feat_data=as.vector(EA_data_import[,1])
EA_data=EA_data_import[,2:length(colnames(EA_data_import))]

#Get only cell lines of interest
EA_data=EA_data[,lines]

#Set rownames for EA data to gene symbol. Should be unique now
rownames(EA_data)=EA_feat_data

#Retrieve corresponding list name
list_name=gene_interest_import[EA_feat_data,2]

#Write results to file
EA_outdata=cbind(rownames(EA_data),EA_data)
write.table(EA_outdata, file="AllGenes_exonarray.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#Get only genes of interest
EA_data=EA_data[EA_feat_data %in% genes_interest,]

#Retrieve corresponding list name
list_name=gene_interest_import[rownames(EA_data),2]
EA_outdata2=cbind(list_name,rownames(EA_data),EA_data)
write.table(EA_outdata2, file="Glycoproteins_exonarray.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



### RNAseq data ###
#1. Gene
RNA_datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Matrix_GeneExpression_v53.txt"
RNA_expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Expressed_GeneExpression_v53.txt"

#Import data
RNA_data_import=read.table(RNA_datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
RNA_exp_status_import=read.table(RNA_expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Break data into features info and expression values
RNA_feat_data=RNA_data_import[,1:4]
RNA_data=RNA_data_import[,5:62]
RNA_exp_status=RNA_exp_status_import[,5:62]


#Get only cell lines of interest
RNA_data=RNA_data[,lines[2:10]] #First line MDAMB468 not in RNAseq data
RNA_exp_status=RNA_exp_status[,lines[2:10]]
RNA_data=log2(RNA_data+1) #Log2 data

#Get only genes of interest
RNA_data=RNA_data[RNA_feat_data[,3] %in% genes_interest,]
RNA_exp_status=RNA_exp_status[RNA_feat_data[,3] %in% genes_interest,]
RNA_feat_data=RNA_feat_data[RNA_feat_data[,3] %in% genes_interest,]

rownames(RNA_data)=RNA_feat_data[,3]
rownames(RNA_exp_status)=RNA_feat_data[,3]

list_name2=gene_interest_import[rownames(RNA_data),2]
RNA_outdata=cbind(list_name2,rownames(RNA_data),RNA_data)
RNA_exp_outdata=cbind(list_name2,rownames(RNA_exp_status),RNA_exp_status)

write.table(RNA_outdata, file="Glycoproteins_RNAseq.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(RNA_exp_outdata, file="Glycoproteins_RNAseq_expressed.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




