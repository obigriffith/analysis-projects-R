dataf <- read.table("erin.rdat")
dataf[,1] <- as.character(dataf[,1])
library(ctest)
get.pval.ttest <- function(dataf,index1,index2,datafilter=as.numeric) {
   f <- function(i) {
     return(t.test(datafilter(dataf[i,index1]),datafilter(dataf[i,index2]))$p.value)
   }
   return(sapply(1:length(dataf[,1]),f))
}
pVal.ttest <- get.pval.ttest(dataf,2:4, 5:7)
print (cbind(dataf,pVal.ttest))
orders <- order(pVal.ttest)
ordered.data <- cbind(dataf[orders,],pVal.ttest[orders])
write(t(as.matrix(ordered.data)),ncolumns=length(dataf)+1,file="erin.out")
