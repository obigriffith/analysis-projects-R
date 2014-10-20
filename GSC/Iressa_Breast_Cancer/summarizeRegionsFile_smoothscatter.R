library(geneplotter)
library(RColorBrewer)

#Set working dir
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/SUM149")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/A431")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/BT474")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/HS578T")

#Specify data file
datafile="SUM149_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha001_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha00001_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha099999_nomin_ALL_fp3_3_1_11.regions"
datafile="SUM149_compare_mode_3_standard_alpha099999_nomin_ALL_fp3_3_2_2.regions"
datafile="SUM149_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_2_2.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_2_2.regions"
datafile="SUM149_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_3_1.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"
datafile="SUM149_compare_mode_3_standard_nofilt_nomin_ALL_fp4_0_6.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp4_0_6.regions"
datafile="SUM149_compare_mode_3_standard_alpha01_nomin_ALL_fp4_0_6.regions"

datafile="A431_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions";
datafile="A431_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_3_1.regions"
datafile="A431_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"

datafile="BT474_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions";
datafile="BT474_compare_mode_3_standard_alpha0999_nomin_ALL_fp3_3_1_11.regions"
datafile="BT474_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"
datafile="BT474_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_3_1.regions"

datafile="HS578T_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_2_2.regions";
datafile="HS578T_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_2_2.regions";
datafile="HS578T_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_3_0.regions";
datafile="HS578T_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_0.regions";
datafile="HS578T_compare_mode_3_standard_nofilt_nomin_ALL_fp3_3_3_1.regions";
datafile="HS578T_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions";


#Specify output file for pdf
#SUM149
pdf("SUM149_diff_peaks_chr1_alpha05_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_alpha001_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_alpha00001_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_alpha099999_nomin_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_alpha099999_nomin_new_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_nofilt_nomin_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr1_nofilt_nomin_smoothscatter_new.pdf")
pdf("SUM149_diff_peaks_chr1_alpha05_nomin_smoothscatter_new.pdf")
pdf("SUM149_diff_peaks_chr1_alpha01_nomin_smoothscatter_new.pdf")

pdf("SUM149_diff_peaks_chr2_alpha05_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr2_alpha001_smoothscatter.pdf")
pdf("SUM149_diff_peaks_chr2_alpha00001_smoothscatter.pdf")

#A431
pdf("A431_diff_peaks_chr1_alpha05_smoothscatter.pdf")
pdf("A431_diff_peaks_chr1_nofilt_nomin_smoothscatter.pdf")

pdf("A431_diff_peaks_chr2_alpha05_smoothscatter.pdf")

#BT474
pdf("BT474_diff_peaks_chr1_alpha05_smoothscatter.pdf")
pdf("BT474_diff_peaks_chr1_alpha05_nomin_smoothscatter.pdf")
pdf("BT474_diff_peaks_chr1_alpha0999_nomin_smoothscatter.pdf")
pdf("BT474_diff_peaks_chr1_nofilt_nomin_smoothscatter.pdf")

#HS578T
pdf("HS578T_diff_peaks_chr1_alpha05_nomin_smoothscatter.pdf")
pdf("HS578T_diff_peaks_chr1_nofilt_nomin_smoothscatter.pdf")
pdf("HS578T_diff_peaks_chr1_alpha05_nomin_smoothscatter_new.pdf")
pdf("HS578T_diff_peaks_chr1_nofilt_nomin_smoothscatter_new.pdf")

#Read in data to R
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Use a single chr as representative example
chr_data=data[data[,1]=="1",]
#chr_data=data[data[,1]=="2",]

#Remove cases with read count in one library but 0 in the other (use this for calculating regression line)
chr_data_sub=chr_data
chr_data_sub=chr_data_sub[chr_data_sub[,5]>0,]
chr_data_sub=chr_data_sub[chr_data_sub[,6]>0,]

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

max=max(chr_data[,5:6])

#Create smooth scatter from all data
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr1, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr1, alpha=0.001", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr1, alpha=0.00001", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr1, alpha=0.99999", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="HS0679 (s.height)", ylab="HS0678 (c.height)", main="Significantly different peaks for SUM149, chr1, no filter", colramp=colors, nbin=275, xlim=c(0,70), ylim=c(0,70)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="HS0679 (s.height)", ylab="HS0678 (c.height)", main="Significantly different peaks for SUM149, chr1, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,70), ylim=c(0,70)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="HS0679 (s.height)", ylab="HS0678 (c.height)", main="Significantly different peaks for SUM149, chr1, alpha=0.01", colramp=colors, nbin=275, xlim=c(0,70), ylim=c(0,70)) 

smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr2, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr2, alpha=0.001", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for SUM149, chr2, alpha=0.00001", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 

smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for A431, chr1, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,33), ylim=c(0,33)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for A431, chr1, nofilter", colramp=colors, nbin=275, xlim=c(0,33), ylim=c(0,33)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for A431, chr2, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 

smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for BT474, chr1, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,max), ylim=c(0,max)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for BT474, chr1, alpha=0.99999", colramp=colors, nbin=275, xlim=c(0,24), ylim=c(0,24)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for BT474, chr1, no filter", colramp=colors, nbin=275, xlim=c(0,21), ylim=c(0,21)) 

smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for HS578T, chr1, alpha=0.05", colramp=colors, nbin=275, xlim=c(0,18), ylim=c(0,18)) 
smoothScatter(x=chr_data[,5], y=chr_data[,6], xlab="s.height", ylab="c.height", main="Significantly different peaks for HS578T, chr1, no filter", colramp=colors, nbin=275, xlim=c(0,18), ylim=c(0,18)) 


#Plot line of best fit (least square regression) through scatterplot
#x=chr_data[,5]
#y=chr_data[,6]
x=chr_data_sub[,5]
y=chr_data_sub[,6]
lm_line=lm(y~x)
abline(lm_line,col="red",lwd=2)

#Add line determined by findPeaks perpendicular least squares fit (e.g., y = (0.63491)x + 0.04483)
abline(a=0.03238, b=0.61893, col="green", lwd=2) #SUM149 chr1
abline(a=0.02540, b=0.56108, col="green", lwd=2) #SUM149 chr1 - Update
abline(a=-0.03515, b=1.15020, col="green", lwd=2) #A431 chr1
abline(a=0.04483, b=0.63491, col="green", lwd=2) #HS578T chr1
abline(a=0.03312, b=0.85326, col="green", lwd=2) #BT474 chr1

#Create confidence intervals. The confidence intervals calculated are extremely small. 
#The prediction (tolerance) intervals seem more reasonable. Not sure if this correct to use though.
#predict_points = seq(1, 24, 1)
predict_points = c(0,sort(unique(x))) #Use this unless there are extreme outliers
#pred.w.clim = predict(lm(y~x), data.frame(x=predict_points), interval="confidence", level = 0.999)
pred.w.clim = predict(lm(y~x), data.frame(x=predict_points), interval="prediction", level = 0.95)
lines(x=predict_points,y=pred.w.clim[,2],lty=2,lwd=2,col="red")
lines(x=predict_points,y=pred.w.clim[,3],lty=2,lwd=2,col="red")

dev.off()

