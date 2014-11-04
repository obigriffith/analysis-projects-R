dataf <- read.csv("/home/obig/clustering/R_analysis/scripts/PearsonMatrix.test.csv", header=F)
dataf[,1] <- as.character(dataf[,1])
ncol <- length(dataf[1,])
library(ctest)

out <- file("/home/obig/clustering/R_analysis/scripts/PearsonMatrix.test.out", "w")
cat("geneA", "geneB", "pearson", "pvalue", "\n", file=out, sep="\t")
for (g1 in 1:length(dataf[,1])){
  for (g2 in 1:length(dataf[,1])) {
    pval <- cor.test(as.numeric(dataf[g1,2:ncol]),as.numeric(dataf[g2,2:ncol]))$p.value
    pearson <- cor(as.numeric(dataf[g1,2:ncol]),as.numeric(dataf[g2,2:ncol]))
    cat(dataf[g1,1], dataf[g2,1], pearson, pval, "\n", file=out, sep="\t")
  }
}
close(out)
