#First, specify to use hexbin library
library("hexbin")

#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/R/R_analysis/sage_freq_vs_ratios_density_power.ps", pointsize=1)
#X11()

#Specify file for data to be read from.
data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_freqs_vs_SAGE_ratios.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_freqs_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_AC_vs_micro_w_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_micro_results_values_only.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_micro_results_values_only_05cutoff.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/bin/Cluster/cluster-1.23/sageresults/final_SAGE_results_values_only.txt",list(x=0))
#data2 <- scan("/home/obig/bin/Cluster/cluster-1.23/results/final_microarray_results_values_only.txt",list(y=0))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x0 <- data1$x
x1 <- data1$y

#Create a hexbin object for the data
hbin <- hexbin(x0,x1,xbins=100)

#Plot the hexbin density plot
plot.hexbin(hbin, xlab="Frequencies",ylab="Ratios", xlim=range(-1,1), ylim=range(-1,1),
      main="SAGE ratios vs frequencies for calculating pearson distances for all pairwise comparisons of 7822 genes", style = "nested.centroids",
      legend = 1, cex=2.5, lcex = 2.5, cex.main = 2.5, cex.lab=2.5, cex.sub=2.5, cex.axis=2.5,
      minarea = 0.05, maxarea = 0.8, trans = NULL, inv = NULL,
      border = FALSE, density = -1,
      colramp = function(n){ LinGray(n,beg = 90,end = 15) },
      verbose = getOption("verbose"))

r <- cor(x0, x1, method=c("pearson"))
print (r)

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