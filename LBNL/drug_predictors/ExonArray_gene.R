##see http://www.aroma-project.org/node/37

#source("http://aroma-project.org/hbLite.R");
# OR source("http://www.braju.com/R/hbLite.R");
#hbInstall("aroma.affymetrix");

setwd('~adaemen/data')

library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

##setup the CDF: TCGA
chipType <- "HuEx-1_0-st-v2"
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="DCCg,Spring2008")
print(cdf)
##Next we setup the CEL set with the above custom CDF:
cs <- AffymetrixCelSet$byName("BCGC_2006", cdf=cdf)
##Need to remove a few bad arrays
bad <- c(grep("BD Control|MDAMB435|HCC1007|HCC1500|DU4475",getNames(cs)))
cs$files<-cs$files[-bad]
print(cs)
##In order to do RMA background correction, we setup a correction method and runs it by:
bc <- RmaBackgroundCorrection(cs, tag="tcga")
csBC <- process(bc,verbose=verbose)
##We then setup a quantile normalization method:
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
print(qn)
csN <- process(qn, verbose=verbose)
## summarization
getCdf(csN)

##To fit a summary of the entire gene (i.e. estimate the overall expression for the gene), do:
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
print(plmTr)
##To fit the PLM to all of the data, do:
fit(plmTr, verbose=verbose)

##To extract the estimates (transcript or probeset) use either extractMatrix() or extractDataFrame() on the ChipEffectSet that corresponds to the plm object
cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr, addNames=TRUE)

##bit of prettification
cn <- colnames(trFit)[-(1:5)]
expr <- log2(extractMatrix(cesTr))
colnames(expr) <- cn
rownames(expr) <- trFit[,1]

write.csv(expr,file='breastExon_genelevel_noHCC1500.csv',quote=F)