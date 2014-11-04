#Load the appropriate libraries
library(affy)
library(gcrma)

#Set working directory
setwd("/projects/02/coexpression/raw_data/h_sapiens/AFFY/GSE2109_CELfiles")

############################just.gcrma method##############################
#Get file list
celfiles=list.files(path="/projects/02/coexpression/raw_data/h_sapiens/AFFY/GSE2109_CELfiles")

#Do gcrma normalization
gcrma_exprset=just.gcrma(filenames=celfiles,normalize=TRUE,type="fullmodel",verbose=TRUE,fast=FALSE,optimize.by="memory")

#extract just the values
gcrma_values=exprs(gcrma_exprset)

#Change directory for writing results
setwd("/home/obig/Projects/sub_space_clustering/expO/normalized")

#write.table(gcrma_values, file = "GSE2109_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
write.table(gcrma_values, file = "GSE2109_1026exps_16Aug2006_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")

######################################################


#####################Annotate Probes to Genes or Proteins################################
#To annotate genes
library("annotate")
library("biomaRt")

#Get latest annotation file for chip of interest
source("http://www.bioconductor.org/biocLite.R")
biocLite("hgu133plus2")
library("hgu133plus2")

#Get probe ids from exprset and remove control probes
probe_ids=geneNames(gcrma_exprset)
probe_ids_nc=probe_ids[1:54613] #Remove controls if desired.

#Create a mart object from ensembl Biomart
mart <- useMart("ensembl", "hsapiens_gene_ensembl")
attrbuts=listAttributes(mart)

########################################Uniprot mapping#########################################
#Get uniprot_swissprot_accessions for all probes
annotations=getBM(attributes=c("affy_hg_u133_plus_2","uniprot_swissprot_accession"), filter="affy_hg_u133_plus_2", values=probe_ids_nc, mart=mart, na.value="NA")   

#Remove failed mappings (i.e. no uniprot found)
annotations_notnull=annotations[which(annotations$uniprot_swissprot_accession!=""),]

#Remove redundant entries (where same probe is mapped to same uniprot - these are pointless)
annotations_notnull_unique=unique(annotations_notnull)

#Get number of mappings for each probe (we want only probes that map to one protein)
probe_counts_table=table(annotations_notnull_unique[,1]) 
probe_counts=probe_counts_table[annotations_notnull_unique[,1]]

#Add number of mappings per probe to annotations data
annotations_notnull_unique_mapcount=cbind(annotations_notnull_unique,probe_counts)

#If we order by probe_id we should see something sensible (i.e. that probes with two mappings to two uniprots have a mapcount of two)
#Check against live BioMart to make sure things look right
annotations_notnull_unique_mapcount[order(annotations_notnull_unique_mapcount[,1]),]

#Remove entries where one probe maps to multiple uniprot ids (still allows many probes to map to one uniprot)
annotations_notnull_unique_unambig=annotations_notnull_unique_mapcount[annotations_notnull_unique_mapcount[,3]==1,][,1:2]

#Write these mappings to a separate file.
write.table(annotations_notnull_unique_unambig, file = "hgu133plus2_Biomart_uniprot_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")

#Use this list of mappings to create an expression matrix with gene/protein ids
mapped=merge(annotations_notnull_unique_unambig,gcrma_values, by.x=1, by.y=0)[2:(length(colnames(gcrma_values))+2)]

#Reduce number of decimals so that file is not so big
mapped_formatted=format(mapped, digits=5)

#Write mapped data to file.
#write.table(mapped_formatted, file = "GSE2109_gcrma_mapped.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
write.table(mapped_formatted, file = "GSE2109_1026exps_16Aug2006_gcrma_uniprot_mapped.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")


########################################Ensembl gene mapping########################################
#Get ensembl_gene_id for all probes
annotations_ENSG=getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"), filter="affy_hg_u133_plus_2", values=probe_ids_nc, mart=mart, na.value="NA")   

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
write.table(annotations_ENSG_notnull_unique_unambig, file = "hgu133plus2_Biomart_ENSG_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")

#Use this list of mappings to create an expression matrix with gene ids
mapped_ENSG=merge(annotations_ENSG_notnull_unique_unambig,gcrma_values, by.x=1, by.y=0)[2:(length(colnames(gcrma_values))+2)]

#Reduce number of decimals so that file is not so big
mapped_ENSG_formatted=format(mapped_ENSG, digits=5)

#Write mapped data to file.
write.table(mapped_ENSG_formatted, file = "GSE2109_1026exps_16Aug2006_gcrma_ENSG_mapped.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "-999", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")



### Unused code####

#############################gcrma method#############################
#Has memory issues (try just.gcrma instead)

#Read in the data
#raw.data=ReadAffy(verbose = TRUE)

##GCRMA normalization
#data.gcrma.norm<-gcrma(raw.data)

##Get the important stuff out of the data - the expression estimates for each array
#gcrma<-exprs(data.gcrma.norm)


#############Function for MvsA plots##########################
mva.plot.fun<-function(v1,v2, filename="plots.pdf"){
  #jpeg(filename, width=960, height = 480) #Can use jpg if you prefer
  pdf(filename)
  M<-v2-v1
  A<-0.5*(v1+v2)
  plot(A,M, xlim=c(-1,15), ylim=c(-15,5))
  abline(-1,0)
  abline(1,0)
  dev.off()
}

mva.plot.fun(gcrma[,1],gcrma[,2], "mva_gcrma.pdf")

