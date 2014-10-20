setwd("C:/Documents and Settings/obig/My Documents/Projects/Gambogic_acid")
datafile="WST1_MIP_v_MIP5FU_w_5FU_GA_allreps_norm_Rdata.txt"
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#MIP5FU tests
#5FU vs DMSO/PBS
MIP5FU_wilcox_5FU_vs_DMSO=wilcox.test(x=data[,5], y=data[,6], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_5FU_vs_DMSO

#GA vs DMSO/PBS
MIP5FU_wilcox_GA_vs_DMSO=wilcox.test(x=data[,5], y=data[,7], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_GA_vs_DMSO

#5FU+GA vs DMSO/PBS
MIP5FU_wilcox_combo_vs_DMSO=wilcox.test(x=data[,5], y=data[,8], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_combo_vs_DMSO

#GA vs 5FU
MIP5FU_wilcox_GA_vs_5FU=wilcox.test(x=data[,6], y=data[,7], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_GA_vs_5FU

#5FU+GA vs 5FU
MIP5FU_wilcox_5FU_vs_combo=wilcox.test(x=data[,6], y=data[,8], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_5FU_vs_combo

#5FU+GA vs GA
MIP5FU_wilcox_GA_vs_combo=wilcox.test(x=data[,7], y=data[,8], alternative="two.sided", paired=TRUE)
MIP5FU_wilcox_GA_vs_combo



#MIP101 tests
#5FU vs DMSO/PBS
MIP101_wilcox_5FU_vs_DMSO=wilcox.test(x=data[,1], y=data[,2], alternative="two.sided", paired=TRUE)
MIP101_wilcox_5FU_vs_DMSO

#GA vs DMSO/PBS
MIP101_wilcox_GA_vs_DMSO=wilcox.test(x=data[,1], y=data[,3], alternative="two.sided", paired=TRUE)
MIP101_wilcox_GA_vs_DMSO

#5FU+GA vs DMSO/PBS
MIP101_wilcox_combo_vs_DMSO=wilcox.test(x=data[,1], y=data[,4], alternative="two.sided", paired=TRUE)
MIP101_wilcox_combo_vs_DMSO

#GA vs 5FU
MIP101_wilcox_GA_vs_5FU=wilcox.test(x=data[,2], y=data[,3], alternative="two.sided", paired=TRUE)
MIP101_wilcox_GA_vs_5FU

#5FU+GA vs 5FU
MIP101_wilcox_5FU_vs_combo=wilcox.test(x=data[,2], y=data[,4], alternative="two.sided", paired=TRUE)
MIP101_wilcox_5FU_vs_combo

#5FU+GA vs GA
MIP101_wilcox_GA_vs_combo=wilcox.test(x=data[,3], y=data[,4], alternative="two.sided", paired=TRUE)
MIP101_wilcox_GA_vs_combo

