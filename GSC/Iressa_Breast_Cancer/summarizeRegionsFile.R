setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/SUM149")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/A431")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/BT474")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Iressa_Breast_Cancer/solexa/findPeaks/HS578T")

datafile="SUM149_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha001_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha00001_ALL_fp3_3_1_11.regions";
datafile="SUM149_compare_mode_3_standard_alpha05_sub_ALL_fp3_3_1_11.regions"
datafile="SUM149_compare_mode_3_standard_alpha00001_sub_ALL_fp3_3_1_11.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_1_11.regions"
datafile="SUM149_compare_mode_3_standard_alpha0999_nomin_ALL_fp3_3_1_11.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"
datafile="SUM149_compare_mode_3_standard_alpha05_nomin_ALL_fp4_0_3.regions"

datafile="A431_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions"
datafile="A431_compare_mode_3_standard_alpha001_ALL_fp3_3_1_11.regions"
datafile="A431_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"

datafile="BT474_compare_mode_3_standard_alpha05_ALL_fp3_3_1_11.regions"
datafile="BT474_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"

datafile="HS578T_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_2_2.regions"
datafile="HS578T_compare_mode_3_standard_alpha05_nomin_ALL_fp3_3_3_1.regions"


data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
chr=unique(data[,1])
chr_ordered=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","MT")
chr_ordered1=c(1,2,3,4,5,6,7,8,9,10,11,12)
chr_ordered2=c(13,14,15,16,17,18,19,20,21,22,"X","Y")
chr_ordered2B=c(13,14,15,16,17,18,19,20,21,22,"X")

pdf ("SUM149_diff_peaks_chr1_12.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha001.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha00001.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha05_subpeaks.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha00001_subpeaks.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha05_nomin.pdf")
pdf ("SUM149_diff_peaks_chr1_12_alpha0999_nomin.pdf")

pdf ("A431_diff_peaks_chr1_12.pdf")
pdf ("A431_diff_peaks_chr1_12_alpha05.pdf")
pdf ("A431_diff_peaks_chr1_12_alpha001.pdf")

pdf ("BT474_diff_peaks_chr1_12_alpha05.pdf")

pdf ("HS578T_diff_peaks_chr1_12_alpha05.pdf")

#Library names
libraries=c("HS0678","HS0679") #SUM149
libraries=c("HS0812","HS0814") #A431
libraries=c("HS0867","HS0869") #HS578T
libraries=c("HS0885","HS0887") #BT474
libraries=c("HS1012","HS1014") #PC9

par(mfrow=c(3,4),oma=c(2,0,2,0))
for (chr in chr_ordered1){
  data_sub=data[data[,1]==chr,]
  max=max(data_sub[,5:6])
  plot(x=data_sub[,5],y=data_sub[,6],type="p", xlab=libraries[2], ylab=libraries[1], main=chr, xlim=c(0,max),ylim=c(0,max))
}
title("FindPeaks Peak Height Comparison", outer=TRUE)
mtext(side=1, datafile, outer=TRUE)
dev.off()

pdf ("SUM149_diff_peaks_chr13_Y.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha001.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha_00001.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha_05_subpeaks.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha_00001_subpeaks.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha_05_nomin.pdf")
pdf ("SUM149_diff_peaks_chr13_Y_alpha_0999_nomin.pdf")

pdf ("A431_diff_peaks_chr13_Y.pdf")
pdf ("A431_diff_peaks_chr13_Y_alpha05.pdf")
pdf ("A431_diff_peaks_chr13_Y_alpha001.pdf")

pdf ("BT474_diff_peaks_chr13_Y_alpha05.pdf")

pdf ("HS578T_diff_peaks_chr13_Y_alpha05.pdf")


par(mfrow=c(3,4),oma=c(0,0,2,0))
for (chr in chr_ordered2){
#for (chr in chr_ordered2B){
  data_sub=data[data[,1]==chr,]
  max=max(data_sub[,5:6])
  plot(x=data_sub[,5],y=data_sub[,6],type="p", xlab=libraries[2], ylab=libraries[1], main=chr, xlim=c(0,max),ylim=c(0,max))
}
title("FindPeaks Peak Height Comparison", outer=TRUE)
mtext(side=1, datafile, outer=TRUE)
dev.off()

