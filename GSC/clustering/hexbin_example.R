#First, specify to use hexbin library
library("hexbin")

#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/tmp/somefile.ps", pointsize=1)

#get data from a file (in this case, a tab-delimited file with two columns of numbers)
data1 <- scan("datafile.txt",list(x=0,y=1))

#Assign data from file to x1 and x2 from data1 file using the x and y labels specified in the list function above
x1 <- data1$x
x2 <- data1$y

#Create a hexbin object for the data
hbin <- hexbin(x1,x2,xbins=100)

#Plot the hexbin density plot
plot.hexbin(hbin, xlab="the x variable",ylab="the y variable", xlim=range(-1,1), ylim=range(-1,1),
      main="this is a graph of x versus y", style = "nested.centroids",
      legend = 1, cex=2.5, lcex = 2.5, cex.main = 2.5, cex.lab=2.5, cex.sub=2.5, cex.axis=2.5,
      minarea = 0.05, maxarea = 0.8, trans = NULL, inv = NULL,
      border = FALSE, density = -1,
      colramp = function(n){ LinGray(n,beg = 90,end = 15) },
      verbose = getOption("verbose"))

##Other commands that might be useful
#hexagons(hbin, colramp=terrain.colors, style="nested.lattice")
#hexagons(hbin, colramp=BTY, maxcnt=100)
#plot.hexbin(hbin, style="nested.lattice", colramp=terrain.colors, legend=1)
#range to use in plot
#rx0 <- range(-1,1)
#rx1 <- range(-1,1)