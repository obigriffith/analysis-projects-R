#First, specify to use hexbin library
library("hexbin")

#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/clustering/R_analysis/misc/Microarray_dataOverlap_vs_Pearson.ps", pointsize=1)
#X11()

#Specify file for data to be read from.
#Internal comparisons
data1 <- scan("/home/sage/clustering/matrix/affy_micro_sage_comp/results/final_micro_results_common2all_w_overlaps_nulls_removed.txt",list(a=0,b=1,y=2,x=3))

#Assign data from files to x1 and y1 from data1 file using the x and y labels specified in the list function above
x1 <- data1$x
y1 <- data1$y


#Create a hexbin object for the data
hbin <- hexbin(x1,y1,xbins=100)

#Plot the hexbin density plot
plot.hexbin(hbin, xlab="# of Overlapping datapoints",ylab="Pearson Correlation Coefficient", 
      main="Number of Overlapping Datapoints vs Pearson Correlation for cDNA microarray", style = "nested.centroids",
      legend = 1, cex=2.5, lcex = 2.5, cex.main = 2.5, cex.lab=2.5, cex.sub=2.5, cex.axis=2.5,
      minarea = 0.05, maxarea = 0.8, trans = NULL, inv = NULL,
      border = FALSE, density = -1,
      colramp = function(n){ LinGray(n,beg = 90,end = 15) },
      verbose = getOption("verbose"))

r <- cor(x1, y1, method=c("pearson"))
p <- cor.test(x1, y1, method = "pearson", alternative = "greater")
print (r)
print (p)

##########################################################
#old commands
##########################################################
#plot(x0,x1, xlab="SAGE", ylab="Microarray", main="SAGE vs Microarray pearson distances for all pairwise comparisons of ~7800 genes")
#plot(rx0,rx1,type="n",xlab="sage",ylab="microarray",main="SAGE vs Microarray pearson distances for all pairwise comparisons of ~7800 genes")
#hexagons(hbin, colramp=terrain.colors, style="nested.lattice")
#hexagons(hbin, colramp=BTY, maxcnt=100)
#plot.hexbin(hbin, style="nested.lattice", colramp=terrain.colors, legend=1)
#range to use in plot
#rx0 <- range(-1,1)
#rx1 <- range(-1,1)

#xlim=range(0,1202), ylim=range(-1,1),