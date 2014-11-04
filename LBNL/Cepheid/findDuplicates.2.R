#Load the appropriate libraries
library(affy)
library(gcrma)
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)
library(genefilter)

#Import all CEL files and look for duplicates by correlation analysis
setwd("/csb/home/obig/Projects/cepheid/processed/customCDF")
duplicates_pdf = "ALL_unfiltered_gcrma_sample_correlations.pdf"
duplicates_zoomed_pdf = "ALL_unfiltered_gcrma_sample_correlations_zoomed.pdf"
duplicates_more_zoomed_pdf = "ALL_unfiltered_gcrma_sample_correlations_more_zoomed.pdf"
duplicates_gt097_pdf = "ALL_unfiltered_gcrma_sample_correlations_gt097.pdf"
duplicates_file = "ALL_unfiltered_gcrma_duplicate_samples.txt"
singles_file = "ALL_unfiltered_gcrma_singles_samples.txt"

#Set CDF to use:
#Use customCDF since number of probes is less than standard
cdf="HGU133A_HS_ENTREZG"

#Read in the data, #This will require a big memory machine
raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="/scratch/obig/data/cepheid/unfiltered/all", cdfname=cdf)

##GCRMA normalization
data.gcrma.norm.ALL=gcrma(raw.data.ALL)

##Get the important stuff out of the data - the expression estimates for each array
gcrma.ALL=exprs(data.gcrma.norm.ALL)

#Preliminary gene filtering - filter out invariant genes. This will make no difference to exact duplicates 
#But, it might help separate re-hybs vs samples that are just really similar. Especially if noise is driving correlation
X=gcrma.ALL
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=gcrma.ALL[filt,] 

#Use this data to look for duplicates in the data
x=as.matrix(filt_Data)
cors=cor(x, method="pearson")
cors_tri=cors
cors_tri[lower.tri(cors_tri, diag=TRUE)]=NA #Set one triangle, plus the diagonal to NA

cors_nodiag=cors
diag(cors_nodiag)=NA

pdf(file=duplicates_pdf)
hist(cors_tri, col="blue", breaks=1000, xlim=c(-0.2,1), main="Sample correlations for all 2174 x 2174 samples, all probes", xlab="Pearson cor")
dev.off()

pdf(file=duplicates_zoomed_pdf)
hist(cors_tri, col="blue", breaks=1000, xlim=c(0.7,1), ylim=c(0,2000), main="Sample correlations for all 2174 x 2174 samples, all probes", xlab="Pearson cor")
dev.off()

pdf(file=duplicates_more_zoomed_pdf)
hist(cors_tri, col="blue", breaks=100, xlim=c(0.8,1), ylim=c(0,2000), main="Sample correlations for all 2174 x 2174 samples, all probes", xlab="Pearson cor")
dev.off()

cors_num=cors_tri[!is.na(cors_tri)]
cors_gt097=cors_num[cors_num>0.828]

pdf(file=duplicates_gt097_pdf)
hist(cors_gt097, col="blue", breaks=100, xlim=c(0.828,1), ylim=c(0,2000), main="Sample correlations for all 2174 x 2174 samples, all probes", xlab="Pearson cor")
dev.off()


#Look at known re-hybed pairs
known_rehybs=matrix(c(
"GSE2990_GSM65752.CEL","GSE7390_GSM178013.CEL",
"GSE2990_GSM65762.CEL","GSE7390_GSM178015.CEL",
"GSE2990_GSM65774.CEL","GSE7390_GSM178026.CEL",
"GSE2990_GSM65775.CEL","GSE7390_GSM178027.CEL",
"GSE2990_GSM65779.CEL","GSE7390_GSM178029.CEL",
"GSE2990_GSM65782.CEL","GSE7390_GSM178030.CEL",
"GSE2990_GSM65786.CEL","GSE7390_GSM178031.CEL",
"GSE2990_GSM65787.CEL","GSE7390_GSM178032.CEL",
"GSE2990_GSM65791.CEL","GSE7390_GSM178033.CEL",
"GSE2990_GSM65792.CEL","GSE7390_GSM178034.CEL",
"GSE2990_GSM65795.CEL","GSE7390_GSM178035.CEL",
"GSE2990_GSM65797.CEL","GSE7390_GSM178039.CEL",
"GSE2990_GSM65798.CEL","GSE7390_GSM178040.CEL",
"GSE2990_GSM65800.CEL","GSE7390_GSM178041.CEL",
"GSE2990_GSM65801.CEL","GSE7390_GSM178042.CEL",
"GSE2990_GSM65805.CEL","GSE7390_GSM178043.CEL",
"GSE2990_GSM65811.CEL","GSE7390_GSM178044.CEL",
"GSE2990_GSM65817.CEL","GSE7390_GSM178045.CEL"
),ncol=2,byrow=TRUE)


known_rehybs_cors=c(
cors["GSE2990_GSM65752.CEL","GSE7390_GSM178013.CEL"],
cors["GSE2990_GSM65762.CEL","GSE7390_GSM178015.CEL"],
cors["GSE2990_GSM65774.CEL","GSE7390_GSM178026.CEL"],
cors["GSE2990_GSM65775.CEL","GSE7390_GSM178027.CEL"],
cors["GSE2990_GSM65779.CEL","GSE7390_GSM178029.CEL"],
cors["GSE2990_GSM65782.CEL","GSE7390_GSM178030.CEL"],
cors["GSE2990_GSM65786.CEL","GSE7390_GSM178031.CEL"],
cors["GSE2990_GSM65787.CEL","GSE7390_GSM178032.CEL"],
cors["GSE2990_GSM65791.CEL","GSE7390_GSM178033.CEL"],
cors["GSE2990_GSM65792.CEL","GSE7390_GSM178034.CEL"],
cors["GSE2990_GSM65795.CEL","GSE7390_GSM178035.CEL"],
cors["GSE2990_GSM65797.CEL","GSE7390_GSM178039.CEL"],
cors["GSE2990_GSM65798.CEL","GSE7390_GSM178040.CEL"],
cors["GSE2990_GSM65800.CEL","GSE7390_GSM178041.CEL"],
cors["GSE2990_GSM65801.CEL","GSE7390_GSM178042.CEL"],
cors["GSE2990_GSM65805.CEL","GSE7390_GSM178043.CEL"],
cors["GSE2990_GSM65811.CEL","GSE7390_GSM178044.CEL"],
cors["GSE2990_GSM65817.CEL","GSE7390_GSM178045.CEL"]
)

#Check how many pairs have high correlation
#use >0.975, Because this captures all but one of the pairs known to be rehybed. One pair was a clear outlier at 0.931
#Or, use >0.99 to just catch the perfect replicates. Can't use 1.0 because sometimes even a perfect correlation appears as 0.999999999
cutoff=0.99
samples=rownames(cors_nodiag)
singles=c()

#First determine number of rows with duplicates to make data matrix
dup_count1=0
for (i in 1:length(samples)){
  dup_test=which(cors_nodiag[i,]>cutoff)
  if (length(dup_test)>0){
    dup_count1=dup_count1+1
  }
}

duplicate_pairs=matrix(data = NA, nrow = dup_count1, ncol = 2)

#Go through and identify sample pairs with cor==1
dup_count2=0
for (i in 1:length(samples)){
  sample=samples[i]
  dup_test=which(cors_nodiag[i,]>cutoff)
  if (length(dup_test)>0){
    dup_count2=dup_count2+1
    dup_samples=paste(names(dup_test),sep="",collapse=" ")
    print(paste(i,sample,dup_samples, sep=" "))
    duplicate_pairs[dup_count2,1]=sample
    duplicate_pairs[dup_count2,2]=dup_samples
  }else{
    singles=c(singles,sample)
  }
}

#Add perfect cors and known rehybs
all_dupes=rbind(duplicate_pairs,known_rehybs)

dup_count2 #Make sure you found all the duplicates
write.table(all_dupes, file=duplicates_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(singles, file=singles_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



