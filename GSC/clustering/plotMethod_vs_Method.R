
#First, specify to use hexbin library
library("hexbin")

#Specify the postscript 'driver' so that output is written to a ps file instead of being displayed by X11 'driver'
postscript("/home/obig/clustering/R_analysis/platform_comps/23sage_28micro_95affy/R_output/randomizations/sample_plots/micro_vs_SAGE_logfreq_noplus1_RANDOM_3way_overlap_density_power.ps", pointsize=1)
#X11()

#Specify file for data to be read from.
#Randomized platform comparison for 23sage_28micro_95affy set
data1 <- scan("/home/obig/clustering/R_analysis/platform_comps/23sage_28micro_95affy/values_only/randomized/micro_vs_SAGE_logfreq_noplus1.random1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/platform_comps/23sage_28micro_95affy/values_only/randomized/AFFY_GEO_07apr04_Ponly_Normal_vs_SAGE_logfreq_noplus1.random1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/platform_comps/23sage_28micro_95affy/values_only/randomized/AFFY_GEO_07apr04_Ponly_Normal_vs_micro.random1",list(x=0,y=1))

#New platform comparisons for 23sage_28micro_95affy set
#data1 <- scan("/home/obig/clustering/R_analysis/23sage_28micro_95affy/values_only/micro_vs_SAGE_logfreq_noplus1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/23sage_28micro_95affy/values_only/AFFY_GEO_07apr04_Ponly_Normal_vs_SAGE_logfreq_noplus1.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/23sage_28micro_95affy/values_only/AFFY_GEO_07apr04_Ponly_Normal_vs_micro.txt",list(x=0,y=1))

#New platform comparisons using updated SAGE and Affy data from GEO
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/AFFY_GEO_07apr04_Ponly_Normal_vs_SAGE_logfreq_noplus1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/AFFY_GEO_07apr04_Ponly_Normal_vs_micro",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/micro_vs_SAGE_logfreq_noplus1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/AFFY_GEO_07apr04_Ponly_vs_SAGE_logfreq_noplus1",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/10sage_100micro_100affy/values_only/AFFY_GEO_07apr04_Ponly_vs_micro",list(x=0,y=1))

#Internal comparisons
#data1 <- scan("/home/sage/clustering/matrix/internal_comparisons/affy/all_GEO/affyA_ln_common2all_cluster.1_vs_2.valuesonly_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/internal_comparisons/microarray/microarray_divide2/microarray.1_vs_microarray.2_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/internal_comparisons/sage/log_freqs_divide2/sage_nofilter_log_freqs.1_vs_log_freqs.2.valuesonly",list(x=0,y=1))

#Sage vs Microarray vs Affy pairwise comparisons for 4170 genes common to all three
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_MPonly_micro_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_MPonly_SAGE_nf_log_freqs_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_Ponly_SAGE_nf_log_freqs_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_Ponly_micro_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/values_only/final_AFFY_GEO_vs_micro_results_common2all.values.nulls_removed",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/values_only/final_AFFY_GEO_vs_SAGE_nf_ratios_results_common2all.values",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/values_only/final_AFFY_GEO_vs_SAGE_nf_log_freqs_results_common2all.values",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_micro_SAGE_nf_log_freqs_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_SAGE_nf_log_freqs_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_micro_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_AFFY_SAGE_nf_ratios_results_common2all_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/clustering/R_analysis/sage_micro_affy/final_micro_SAGE_nf_ratios_results_common2all_nulls_removed.txt",list(x=0,y=1))

##Sage vs Microarray pairwise comparisons for 7822 genes common to both
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_freqs_vs_SAGE_ratios.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/SAGE_micro_nofilter_freqs_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_AC_vs_micro_w_nulls_removed.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_micro_results_values_only.txt",list(x=0,y=1))
#data1 <- scan("/home/sage/clustering/matrix/sage_micro_comp/final_SAGE_micro_results_values_only_05cutoff.txt",list(x=0,y=1))
#data1 <- scan("/home/obig/bin/Cluster/cluster-1.23/sageresults/final_SAGE_results_values_only.txt",list(x=0))
#data2 <- scan("/home/obig/bin/Cluster/cluster-1.23/results/final_microarray_results_values_only.txt",list(y=0))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x1 <- data1$x
x2 <- data1$y


#Create a hexbin object for the data
hbin <- hexbin(x1,x2,xbins=100)

#Plot the hexbin density plot
plot.hexbin(hbin, xlab="micro (28)",ylab="SAGE (GEO, logfreq, noplus1, 23)", xlim=range(-1,1), ylim=range(-1,1),
      main="micro (28 overlap) vs SAGE (23 overlap) randomized data", style = "nested.centroids",
      legend = 1, cex=2.5, lcex = 2.5, cex.main = 2.5, cex.lab=2.5, cex.sub=2.5, cex.axis=2.5,
      minarea = 0.05, maxarea = 0.8, trans = NULL, inv = NULL,
      border = FALSE, density = -1,
      colramp = function(n){ LinGray(n,beg = 90,end = 15) },
      verbose = getOption("verbose"))

r <- cor(x1, x2, method=c("pearson"))
p <- cor.test(x1, x2, method = "pearson", alternative = "greater")
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