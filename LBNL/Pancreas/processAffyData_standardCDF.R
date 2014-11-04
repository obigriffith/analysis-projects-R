#Load the appropriate libraries
library(affy)
library(gcrma)
library(hgu133plus2.db)

#Set working directory
setwd("/csb/home/obig/Projects/pancreas/processed")

#Read in the data
raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/csb/home/obig/Projects/pancreas/Collisson_Hu133PLus2_CEL")

##GCRMA normalization
data.gcrma.norm.ALL=gcrma(raw.data.ALL)

##Get the important stuff out of the data - the expression estimates for each array
gcrma.ALL=exprs(data.gcrma.norm.ALL)

#Remove control probes
gcrma.ALL=gcrma.ALL[1:54613,] #Remove Affy control probes, affy CDF

#Format to 5 decimal places
gcrma.ALL=format(gcrma.ALL, digits=5)

#Map probes to gene symbols
#To see all mappings for Entrez gene db associated with customCDF
ls("package:hgu133plus2.db")
probes.ALL=row.names(gcrma.ALL)
symbol.ALL = unlist(mget(probes.ALL, hgu133plus2SYMBOL))
ID.ALL = unlist(mget(probes.ALL, hgu133plus2ENTREZID))

gcrma.ALL=cbind(probes.ALL,ID.ALL,symbol.ALL,gcrma.ALL)

#Write GCRMA-normalized, mapped data to file
write.table(gcrma.ALL, file = "ALL_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

