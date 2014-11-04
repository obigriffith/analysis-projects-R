#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/clustering/R_analysis/freq_dists/Affy_cum_freq_dist_hist.ps", pointsize=1)


#Specify file for data to be read from.
#new analysis
#data1 <- scan("/home/sage/clustering/matrix/affy_micro_sage_comp/040407/results/for_coexpression_db/affy_only_gt100.txt",list(x=2))
data1 <- scan("/home/obig/clustering/R_analysis/temp/affy_test_data.txt",list(x=2))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x0 <- data1$x

cumhist = function(x)
{
  Z = hist( x0 , plot=F )
  plot(1:length(Z$counts),cumsum(Z$counts)/length(x)*100,type="b",axes=F,
       ylab="%",xlab="")
  axis(1,at=1:length(Z$counts),labels=round(1/Z$mids,digits=0))
  axis(2)
}

#hist(x0, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="freq. dist. of pears. values calcs. for cDNA (28 overlap)", plot=TRUE)

#Other options
#breaks="Sturges" Comes up with some reasonable choice of breaks.