##see http://www.aroma-project.org/node/37

#source("http://aroma-project.org/hbLite.R");
# OR source("http://www.braju.com/R/hbLite.R");
#hbInstall("aroma.affymetrix");

setwd('~adaemen/data')

library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)


##setup the CDF: Affy for transcript and exon summaries
chipType <- "HuEx-1_0-st-v2"
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP")
print(cdf)
##Next we setup the CEL set with the above custom CDF:
cs <- AffymetrixCelSet$byName("BCGC_2006", cdf=cdf)
##Need to remove a few bad arrays
bad <- c(grep("BD Control|MDAMB435|HCC1007|HCC1500|DU4475",getNames(cs)))
cs$files<-cs$files[-bad]
print(cs)
##In order to do RMA background correction, we setup a correction method and runs it by:
bc <- RmaBackgroundCorrection(cs, tag="coreR3")
csBC <- process(bc,verbose=verbose)
##We then setup a quantile normalization method:
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
print(qn)
csN <- process(qn, verbose=verbose)
## summarization
getCdf(csN)


##To fit exon-by-exon, do:
plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
print(plmEx)
fit(plmEx, verbose=verbose)
cesEx <- getChipEffectSet(plmEx)
exFit <- extractDataFrame(cesEx, addNames=TRUE)

cnEx <- colnames(exFit)[-(1:5)]
exprEx <- log2(extractMatrix(cesEx))
colnames(exprEx) <- cnEx
rownames(exprEx) <- exFit[,2]

write.csv(exprEx,file='breastExon_exonlevel_core_noHCC1500.csv',quote=F)
