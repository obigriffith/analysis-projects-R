library("hexbin")
#postscript("sage_vs_micro_density.ps", pointsize=1)
X11()
dataf <- read.table("/home/obig/R/scripts/test.data")

x0 <- as.numeric(dataf[,1])
x1 <- as.numeric(dataf[,2])

rx0 <- range(x0)
rx1 <- range(x1)
hbin <- hexbin(x0,x1,xbins=10)

plot(rx0,rx1,type="n",xlab="sage",ylab="microarray")
hexagons(hbin)



