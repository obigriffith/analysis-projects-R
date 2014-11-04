#use one gene per line, one experiment per column
#edit Ngenes and Nexps accordingly
#use NA for non-available data

Ngenes <- 4
Nexps <- 7

data <- matrix(scan("/home/obig/clustering/R_analysis/scripts/erin3.mat", n=Ngenes*Nexps), Ngenes, Nexps, byrow=TRUE)
library(ctest)

out <- file("/home/obig/clustering/R_analysis/scripts/erin3.out", "w")
cat("geneA", "geneB", "pearson", "pvalue", "\n", file=out, sep="\t")
for (g1 in 1:Ngenes){
  for (g2 in 1:Ngenes) {
    if (g2 > g1) {
      pval <- cor.test(data[g1,],data[g2,])$p.value
      pearson <- cor(data[g1,],data[g2,],use="complete.obs")
      cat(g1, g2, pearson, pval, "\n", file=out, sep="\t")
    }
  }
}
close(out)
