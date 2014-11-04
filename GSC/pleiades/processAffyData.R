#Load the appropriate libraries
library(affy)
library(gcrma)
library("annotate")
library("biomaRt")

#Get latest annotation file for platform of interest (e.g. hgu133plus2, mgu74av2, etc)
source("http://www.bioconductor.org/biocLite.R")
biocLite("mgu74av2")
library("mgu74av2")

#Set working directory
setwd("/home/obig/Projects/pleiades/GSE3594/CELfiles")

############################just.gcrma method##############################
#Get file list
celfiles=list.files(path="/home/obig/Projects/pleiades/GSE3594/CELfiles")

#Do gcrma normalization.  Optimize for speed or memory depending on size of dataset
#gcrma_exprset=just.gcrma(filenames=celfiles,normalize=TRUE,type="fullmodel",verbose=TRUE,fast=FALSE,optimize.by="memory")
gcrma_exprset=just.gcrma(filenames=celfiles,normalize=TRUE,type="fullmodel",verbose=TRUE,fast=FALSE,optimize.by="speed")

#extract just the values
gcrma_values=exprs(gcrma_exprset)

#Change directory for writing results
setwd("/home/obig/Projects/pleiades/GSE3594/normalized")

#Reduce number of decimals so that file is not so big
gcrma_values_formatted=format(gcrma_values, digits=5)

write.table(gcrma_values_formatted, file = "GSE3594_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
######################################################


#####################Annotate Probes to Genes or Proteins################################
#Get probe ids from exprset and remove control probes
probe_ids=geneNames(gcrma_exprset)

#Remove controls if desired. Make sure to use correct range for platform in question
probe_ids_nc=probe_ids[1:12422] 

#Create a mart object from ensembl Biomart
mart <- useMart("ensembl", "mmusculus_gene_ensembl")
attrbuts=listAttributes(mart)


########################################Ensembl gene mapping########################################
#Get ensembl_gene_id for all probes
annotations_ENSG=getBM(attributes=c("affy_mg_u74av2","ensembl_gene_id"), filter="affy_mg_u74av2", values=probe_ids_nc, mart=mart, na.value="NA")   

#Remove failed mappings (i.e. no ENSG found)
annotations_ENSG_notnull=annotations_ENSG[which(annotations_ENSG$ensembl_gene_id!=""),]

#Remove redundant entries (where same probe is mapped to same ENSG - these are pointless)
annotations_ENSG_notnull_unique=unique(annotations_ENSG_notnull)

#Get number of mappings for each probe (we want only probes that map to one ENSG)
probe_counts_table_ENSG=table(annotations_ENSG_notnull_unique[,1]) 
probe_counts_ENSG=probe_counts_table_ENSG[annotations_ENSG_notnull_unique[,1]]

#Add number of mappings per probe to annotations data
annotations_ENSG_notnull_unique_mapcount=cbind(annotations_ENSG_notnull_unique,probe_counts_ENSG)

#If we order by probe_id we should see something sensible (i.e. that probes with two mappings to two uniprots have a mapcount of two)
#Check against live BioMart to make sure things look right
annotations_ENSG_notnull_unique_mapcount[order(annotations_ENSG_notnull_unique_mapcount[,1]),]

#Remove entries where one probe maps to multiple ENSG ids (still allows many probes to map to one ENSG)
annotations_ENSG_notnull_unique_unambig=annotations_ENSG_notnull_unique_mapcount[annotations_ENSG_notnull_unique_mapcount[,3]==1,][,1:2]

#Write these mappings to a separate file.
write.table(annotations_ENSG_notnull_unique_unambig, file = "mgu74av2_Biomart_ENSG_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")

#Use this list of mappings to create an expression matrix with gene ids
mapped_ENSG=merge(annotations_ENSG_notnull_unique_unambig,gcrma_values, by.x=1, by.y=0)[2:(length(colnames(gcrma_values))+2)]

#Reduce number of decimals so that file is not so big
mapped_ENSG_formatted=format(mapped_ENSG, digits=5)

#Write mapped data to file.
write.table(mapped_ENSG_formatted, file = "GSE3594_gcrma_ENSG_mapped.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
