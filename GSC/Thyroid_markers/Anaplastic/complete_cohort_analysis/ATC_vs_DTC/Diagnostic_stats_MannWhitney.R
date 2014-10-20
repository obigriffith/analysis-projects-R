library(multtest)

#Set working directory and filenames for Input/output
#setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimfoci_vs_DTC")
#setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimFociNoDiff_vs_DTC")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCdiff_Foci_vs_DTC")

#datafile="ATCprimfoci_vs_DTC.txt"
#datafile="ATCprimFociNoDiff_vs_DTC.txt"
datafile="ATCdiff_Foci_vs_DTC.txt"

#outfile="ATCprimfoci_vs_DTC_ungrouped_MannWhitney.txt"
#outfile="ATCprimFociNoDiff_vs_DTC_ungrouped_MannWhitney.txt"
outfile="ATCdiff_Foci_vs_DTC_ungrouped_MannWhitney.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
path=data[,7]

#Get marker data
marker_data=data[,10:63] #ungrouped values

marker_names=colnames(marker_data)

#Get diagnostic target variable and specify as factor/categorical
target=path
target[target==0]="ATC"
target[target==1]="DTC"
target=as.factor(target)

#Create array to store results for diagnostic contingency table stats
wilcox_results = array(0, dimnames = list(marker_names, c("wilcox_pvalue", "W", "ATC_mean_rank", "DTC_mean_rank", "Direction")), dim=c(length(marker_names),5))

#For each marker perform the Mann-Whitney U test for marker score vs. pathology (ATC vs. DTC)
for (i in 1:length(marker_data)){
   #Summarize data for marker vs diagnosis
   #Calculate mean ranks for ATC vs DTC so that direction of change can be determined
   ranked_marker_data=rank(marker_data[,i], na.last="keep")
   ATC_mean_rank=mean(ranked_marker_data[target=="ATC"], na.rm=TRUE)
   DTC_mean_rank=mean(ranked_marker_data[target=="DTC"], na.rm=TRUE)
   wilcox_results[i,"ATC_mean_rank"]=ATC_mean_rank
   wilcox_results[i,"DTC_mean_rank"]=DTC_mean_rank
   #Determine if marker score has increased or decreased
   if(DTC_mean_rank-ATC_mean_rank>0){
      direction="UP"
   }
   if(DTC_mean_rank-ATC_mean_rank<0){
      direction="DOWN"
   }
   if(DTC_mean_rank-ATC_mean_rank==0){
      direction=NA
   }
   wilcox_results[i,"Direction"]=direction

   #Perform Mann-Whitney U test
   ATC_scores=marker_data[target=="ATC",i]
   DTC_scores=marker_data[target=="DTC",i]
   wilcox_result=wilcox.test(x=ATC_scores, y=DTC_scores, alternative="two.sided", paired=FALSE)
   wilcox_results[i,"wilcox_pvalue"]=wilcox_result$p.value
   wilcox_results[i,"W"]=wilcox_result$statistic
}

#Correct p-values
wilcox_pvalues=as.numeric(wilcox_results[,"wilcox_pvalue"])
wilcox_pvalues_adj=mt.rawp2adjp(wilcox_pvalues, proc=c("Bonferroni","BH"))
wilcox_pvalues_adj_orig_order=wilcox_pvalues_adj$adjp[order(wilcox_pvalues_adj$index),]
wilcox_results=cbind(wilcox_results, wilcox_pvalues_adj_orig_order[,2:3])

#Write complete results to file
write.table (wilcox_results, sep="\t", file=outfile)

