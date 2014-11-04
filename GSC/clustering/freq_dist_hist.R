#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/clustering/R_analysis/23sage_28micro_95affy/R_output/freq_dists/micro_pearson_distribution.ps", pointsize=1)
#postscript("/home/obig/R/R_analysis/test_pears_distribution.ps", pointsize=1)

#Specify file for data to be read from.
#new analysis
data1 <- scan("/home/obig/clustering/R_analysis/23sage_28micro_95affy/values_only/micro_vs_SAGE_logfreq_noplus1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/AFFY_GEO_07apr04_Ponly_Normal_vs_micro",list(x=0,y=1))

#old analysis
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_log_freqs_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_freqs_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_dror_cut05_vs_micro_w_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_AC_vs_micro_w_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_freqs_vs_SAGE_ratios.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/R/scripts/test.data2",list(x=0,y=1))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x0 <- data1$x
y0 <- data1$y

hist(x0, breaks=201,freq=TRUE, col="blue", xlim=range(-1,1), xlab="Pearson correlation coefficient (r)", main="freq. dist. of pears. values calcs. for cDNA (28 overlap)", plot=TRUE)


#Other options
#breaks="Sturges" Comes up with some reasonable choice of breaks.