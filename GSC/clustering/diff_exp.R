library(Biobase, warn.conflicts = FALSE)
library(affy, warn.conflicts=FALSE)
library(affydata)
library(stats)
library(annotate)
library(hgu95av2)

#Read CEL file list into R
spikein <- ReadAffy("FILE1.CEL", "FILE2.CEL")

#Do some Quality Control
#Check distribution of log intensities
pops <- pData(spikein)[, 1] + 2
hist(spikein, col = pops, type = "l")

#Do some pre-processing
#This converts probe-level data into expression measures (whatever that means)
eset <- rma(spikein)

#The following should tell you how many populations and replicates you have 
eset$population

#To get average intensities, average log ratios, t-tests, and p-values using the t.test function
Index1 <- which(eset$population == 0)
Index2 <- which(eset$population == 1)
scores <- esApply(eset, 1, function(x){
tmp <- t.test(x[Index2], x[Index1], var.equal = TRUE)
c(mean(tmp$estimate), -diff(tmp$estimate), tmp$statistic,
tmp$p.value)
})


abatch <- ReadAffy(celfile.path=celpath)