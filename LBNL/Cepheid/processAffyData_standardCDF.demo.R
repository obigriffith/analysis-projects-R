#Load the appropriate libraries
library(affy)
library(gcrma)
library(hgu133a.db) #cdfname="HGU133A"

#Set working directory - where results file will go
setwd("/csb/home/obig/Projects/cepheid/processed_final2/test_survival/standardCDF")

#Set CDF to use:
cdf="HGU133A"

#Read in the data
raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/test_survival") #To use standard CDF

##GCRMA normalization
data.gcrma.norm.ALL=gcrma(raw.data.ALL)

##Get the important stuff out of the data - the expression estimates for each array
gcrma.ALL=exprs(data.gcrma.norm.ALL)

#Remove control probes
gcrma.ALL=gcrma.ALL[1:22215,] #Remove Affy control probes, affy CDF

#Format to 5 decimal places
gcrma.ALL=format(gcrma.ALL, digits=5)

#Map probes to gene symbols
#To see all mappings for Entrez gene db associated with customCDF
ls("package:hgu133a.db") #standrardCDF

probes.ALL=row.names(gcrma.ALL)

#standardDCF
symbol.ALL = unlist(mget(probes.ALL, hgu133aSYMBOL))
ID.ALL = unlist(mget(probes.ALL, hgu133aENTREZID))

gcrma.ALL=cbind(probes.ALL,ID.ALL,symbol.ALL,gcrma.ALL)

#Write GCRMA-normalized, mapped data to file
write.table(gcrma.ALL, file = "ALL_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

