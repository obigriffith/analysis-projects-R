#Load the appropriate libraries
library(affy)
library(gcrma)
#library(hgu133ahsensgcdf) #cdfname="HGU133A_HS_ENSG"
#library(hgu133ahsensgprobe)
#library(hgu133ahsenstcdf) #cdfname="HGU133A_HS_ENST"
#library(hgu133ahsenstprobe) 
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)
#library(hgu133a.db) #cdfname="HGU133A"

#Set working directory
#setwd("/csb/home/obig/Projects/cepheid/processed/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed/standardCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed2/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed2/standardCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final/standardCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final2/train_survival/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final2/train_survival/standardCDF")

setwd("/csb/home/obig/Projects/cepheid/processed_final2/test_survival/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final2/test_survival/standardCDF")

#setwd("/csb/home/obig/Projects/cepheid/processed_final2/test_train_survival/customCDF")
#setwd("/csb/home/obig/Projects/cepheid/processed_final2/test_train_survival/standardCDF")


#Set CDF to use:
#cdf="HGU133A_HS_ENSG"
#cdf="HGU133A_HS_ENST"
cdf="HGU133A_HS_ENTREZG"
#cdf="HGU133A"

#Read in the data
#raw.data.GSE11121=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE11121", cdfname=cdf)
#raw.data.GSE12093=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE12093", cdfname=cdf)
#raw.data.GSE17705=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE17705", cdfname=cdf)
#raw.data.GSE2034=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE2034", cdfname=cdf)
#raw.data.GSE2990=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE2990", cdfname=cdf)
#raw.data.GSE4922=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE4922", cdfname=cdf)
#raw.data.GSE6532=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE6532", cdfname=cdf)
#raw.data.GSE7390=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/GSE7390", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/train_test_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered/train_test_survival") #To use standard CDF
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered2/train_test_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered2/train_test_survival") #To use standard CDF
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final/train_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final/train_survival") #To use standard CDF
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/train_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/train_survival") #To use standard CDF

raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/test_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/test_survival") #To use standard CDF

#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/test_train_survival", cdfname=cdf)
#raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/filtered_final2/test_train_survival") #To use standard CDF

##GCRMA normalization
#data.gcrma.norm.GSE11121=gcrma(raw.data.GSE11121)
#data.gcrma.norm.GSE12093=gcrma(raw.data.GSE12093)
#data.gcrma.norm.GSE17705=gcrma(raw.data.GSE17705)
#data.gcrma.norm.GSE2034=gcrma(raw.data.GSE2034)
#data.gcrma.norm.GSE2990=gcrma(raw.data.GSE2990)
#data.gcrma.norm.GSE4922=gcrma(raw.data.GSE4922)
#data.gcrma.norm.GSE6532=gcrma(raw.data.GSE6532)
#data.gcrma.norm.GSE7390=gcrma(raw.data.GSE7390)
data.gcrma.norm.ALL=gcrma(raw.data.ALL)

##Get the important stuff out of the data - the expression estimates for each array
#gcrma.GSE11121=exprs(data.gcrma.norm.GSE11121)
#gcrma.GSE12093=exprs(data.gcrma.norm.GSE12093)
#gcrma.GSE17705=exprs(data.gcrma.norm.GSE17705)
#gcrma.GSE2034=exprs(data.gcrma.norm.GSE2034)
#gcrma.GSE2990=exprs(data.gcrma.norm.GSE2990)
#gcrma.GSE4922=exprs(data.gcrma.norm.GSE4922)
#gcrma.GSE6532=exprs(data.gcrma.norm.GSE6532)
#gcrma.GSE7390=exprs(data.gcrma.norm.GSE7390)
gcrma.ALL=exprs(data.gcrma.norm.ALL)

#Remove control probes
#gcrma.GSE11121=gcrma.GSE11121[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE12093=gcrma.GSE12093[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE17705=gcrma.GSE17705[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE2034=gcrma.GSE2034[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE2990=gcrma.GSE2990[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE4922=gcrma.GSE4922[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE6532=gcrma.GSE6532[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.GSE7390=gcrma.GSE7390[1:12065,] #Remove Affy control probes, custom CDF
gcrma.ALL=gcrma.ALL[1:12065,] #Remove Affy control probes, custom CDF
#gcrma.ALL=gcrma.ALL[1:22215,] #Remove Affy control probes, affy CDF

#Format to 5 decimal places
#gcrma.GSE11121=format(gcrma.GSE11121, digits=5)
#gcrma.GSE12093=format(gcrma.GSE12093, digits=5)
#gcrma.GSE17705=format(gcrma.GSE17705, digits=5)
#gcrma.GSE2034=format(gcrma.GSE2034, digits=5)
#gcrma.GSE2990=format(gcrma.GSE2990, digits=5)
#gcrma.GSE4922=format(gcrma.GSE4922, digits=5)
#gcrma.GSE6532=format(gcrma.GSE6532, digits=5)
#gcrma.GSE7390=format(gcrma.GSE7390, digits=5)
gcrma.ALL=format(gcrma.ALL, digits=5)

#Map probes to gene symbols
#To see all mappings for Entrez gene db associated with customCDF
ls("package:hgu133ahsentrezg.db") #customCDF
#ls("package:hgu133a.db") #standrardCDF

#probes.GSE11121=row.names(gcrma.GSE11121)
#probes.GSE12093=row.names(gcrma.GSE12093)
#probes.GSE17705=row.names(gcrma.GSE17705)
#probes.GSE2034=row.names(gcrma.GSE2034)
#probes.GSE2990=row.names(gcrma.GSE2990)
#probes.GSE4922=row.names(gcrma.GSE4922)
#probes.GSE6532=row.names(gcrma.GSE6532)
#probes.GSE7390=row.names(gcrma.GSE7390)
probes.ALL=row.names(gcrma.ALL)

#customCDF
#symbol.GSE11121 = unlist(mget(probes.GSE11121, hgu133ahsentrezgSYMBOL))
#symbol.GSE12093 = unlist(mget(probes.GSE12093, hgu133ahsentrezgSYMBOL))
#symbol.GSE17705 = unlist(mget(probes.GSE17705, hgu133ahsentrezgSYMBOL))
#symbol.GSE2034 = unlist(mget(probes.GSE2034, hgu133ahsentrezgSYMBOL))
#symbol.GSE2990 = unlist(mget(probes.GSE2990, hgu133ahsentrezgSYMBOL))
#symbol.GSE4922 = unlist(mget(probes.GSE4922, hgu133ahsentrezgSYMBOL))
#symbol.GSE6532 = unlist(mget(probes.GSE6532, hgu133ahsentrezgSYMBOL))
#symbol.GSE7390 = unlist(mget(probes.GSE7390, hgu133ahsentrezgSYMBOL))
symbol.ALL = unlist(mget(probes.ALL, hgu133ahsentrezgSYMBOL))

#standardDCF
#symbol.GSE11121 = unlist(mget(probes.GSE11121, hgu133aSYMBOL))
#symbol.GSE12093 = unlist(mget(probes.GSE12093, hgu133aSYMBOL))
#symbol.GSE17705 = unlist(mget(probes.GSE17705, hgu133aSYMBOL))
#symbol.GSE2034 = unlist(mget(probes.GSE2034, hgu133aSYMBOL))
#symbol.GSE2990 = unlist(mget(probes.GSE2990, hgu133aSYMBOL))
#symbol.GSE4922 = unlist(mget(probes.GSE4922, hgu133aSYMBOL))
#symbol.GSE6532 = unlist(mget(probes.GSE6532, hgu133aSYMBOL))
#symbol.GSE7390 = unlist(mget(probes.GSE7390, hgu133aSYMBOL))
#symbol.ALL = unlist(mget(probes.ALL, hgu133aSYMBOL))

#customCDF
#ID.GSE11121 = unlist(mget(probes.GSE11121, hgu133ahsentrezgENTREZID))
#ID.GSE12093 = unlist(mget(probes.GSE12093, hgu133ahsentrezgENTREZID))
#ID.GSE17705 = unlist(mget(probes.GSE17705, hgu133ahsentrezgENTREZID))
#ID.GSE2034 = unlist(mget(probes.GSE2034, hgu133ahsentrezgENTREZID))
#ID.GSE2990 = unlist(mget(probes.GSE2990, hgu133ahsentrezgENTREZID))
#ID.GSE4922 = unlist(mget(probes.GSE4922, hgu133ahsentrezgENTREZID))
#ID.GSE6532 = unlist(mget(probes.GSE6532, hgu133ahsentrezgENTREZID))
#ID.GSE7390 = unlist(mget(probes.GSE7390, hgu133ahsentrezgENTREZID))
ID.ALL = unlist(mget(probes.ALL, hgu133ahsentrezgENTREZID))

#standardDCF
#ID.GSE11121 = unlist(mget(probes.GSE11121, hgu133aENTREZID))
#ID.GSE12093 = unlist(mget(probes.GSE12093, hgu133aENTREZID))
#ID.GSE17705 = unlist(mget(probes.GSE17705, hgu133aENTREZID))
#ID.GSE2034 = unlist(mget(probes.GSE2034, hgu133aENTREZID))
#ID.GSE2990 = unlist(mget(probes.GSE2990, hgu133aENTREZID))
#ID.GSE4922 = unlist(mget(probes.GSE4922, hgu133aENTREZID))
#ID.GSE6532 = unlist(mget(probes.GSE6532, hgu133aENTREZID))
#ID.GSE7390 = unlist(mget(probes.GSE7390, hgu133aENTREZID))
#ID.ALL = unlist(mget(probes.ALL, hgu133aENTREZID))

#gcrma.GSE11121=cbind(probes.GSE11121,ID.GSE11121,symbol.GSE11121,gcrma.GSE11121)
#gcrma.GSE12093=cbind(probes.GSE12093,ID.GSE12093,symbol.GSE12093,gcrma.GSE12093)
#gcrma.GSE17705=cbind(probes.GSE17705,ID.GSE17705,symbol.GSE17705,gcrma.GSE17705)
#gcrma.GSE2034=cbind(probes.GSE2034,ID.GSE2034,symbol.GSE2034,gcrma.GSE2034)
#gcrma.GSE2990=cbind(probes.GSE2990,ID.GSE2990,symbol.GSE2990,gcrma.GSE2990)
#gcrma.GSE4922=cbind(probes.GSE4922,ID.GSE4922,symbol.GSE4922,gcrma.GSE4922)
#gcrma.GSE6532=cbind(probes.GSE6532,ID.GSE6532,symbol.GSE6532,gcrma.GSE6532)
#gcrma.GSE7390=cbind(probes.GSE7390,ID.GSE7390,symbol.GSE7390,gcrma.GSE7390)
gcrma.ALL=cbind(probes.ALL,ID.ALL,symbol.ALL,gcrma.ALL)

#Write GCRMA-normalized, mapped data to file
#write.table(gcrma.GSE11121, file = "GSE11121_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE12093, file = "GSE12093_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE17705, file = "GSE17705_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE2034, file = "GSE2034_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE2990, file = "GSE2990_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE4922, file = "GSE4922_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE6532, file = "GSE6532_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(gcrma.GSE7390, file = "GSE7390_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(gcrma.ALL, file = "ALL_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

