library(multtest)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimfoci_vs_DTC")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimFociNoDiff_vs_DTC")
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCdiff_Foci_vs_DTC")

datafile="ATCprimfoci_vs_DTC.txt"
datafile="ATCprimFociNoDiff_vs_DTC.txt"
datafile="ATCdiff_Foci_vs_DTC.txt"

#outfile="ATCprimfoci_vs_DTC_group1_diagtests.txt"
#outfile="ATCprimFociNoDiff_vs_DTC_group1_diagtests.txt"
#outfile="ATCdiff_Foci_vs_DTC_group1_diagtests.txt"
outfile="ATCprimfoci_vs_DTC_group2_diagtests.txt"
outfile="ATCprimFociNoDiff_vs_DTC_group2_diagtests.txt"
outfile="ATCdiff_Foci_vs_DTC_group2_diagtests.txt"


#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

path=data[,7]

#Get marker data
marker_data=data[,64:117] #group 1 values
marker_data=data[,118:171] #group 2 values

marker_names=colnames(marker_data)

#Get diagnostic target variable and specify as factor/categorical
target=path
target[target==0]="ATC"
target[target==1]="DTC"
target=as.factor(target)

#Make sure predictor data is correctly recognized as categorical
for(i in 1:length(marker_data)){
marker_data[,i]=as.factor(marker_data[,i])
}

#Create array to store results for diagnostic contingency table stats
cont_table_results = array(0, dimnames = list(marker_names, c("ATC_Neg", "ATC_Pos", "DTC_Neg", "DTC_Pos", "ATC_Pos_perc", "DTC_Pos_perc", "cont_table_pvalue", "test")), dim=c(length(marker_names),8))

for (i in 1:length(marker_data)){
   #Summarize data for marker vs diagnosis
   ATC_scores=marker_data[target=="ATC",i]
   DTC_scores=marker_data[target=="DTC",i]
   ATC_scores_noNA=ATC_scores[!is.na(ATC_scores)]
   DTC_scores_noNA=DTC_scores[!is.na(DTC_scores)]
   ATC_Neg=length(ATC_scores_noNA[ATC_scores_noNA==0])
   ATC_Pos=length(ATC_scores_noNA[ATC_scores_noNA==1])
   DTC_Neg=length(DTC_scores_noNA[DTC_scores_noNA==0])
   DTC_Pos=length(DTC_scores_noNA[DTC_scores_noNA==1])
   ATC_Pos_perc=(ATC_Pos/(ATC_Pos+ATC_Neg))*100
   DTC_Pos_perc=(DTC_Pos/(DTC_Pos+DTC_Neg))*100
   cont_table_results[i,"ATC_Neg"]=ATC_Neg
   cont_table_results[i,"ATC_Pos"]=ATC_Pos
   cont_table_results[i,"DTC_Neg"]=DTC_Neg
   cont_table_results[i,"DTC_Pos"]=DTC_Pos
   cont_table_results[i,"ATC_Pos_perc"]=ATC_Pos_perc
   cont_table_results[i,"DTC_Pos_perc"]=DTC_Pos_perc

   #Skip any markers which do not have at least two factors (i.e. if constant)
   if (nlevels(marker_data[,i])<2){
      cont_table_results[i,"cont_table_pvalue"]=NA
      cont_table_results[i,"test"]=NA
      next
   }

   #Calculate Pearson Chi-square for all markers vs diagnosis
   chisquare_result=chisq.test(x=marker_data[,i], y=target)
   observed=chisquare_result$observed
   expected=chisquare_result$expected

   #If one or more cells in the table have expected count less than or equal to 5 use a Fisher test instead of Pearson Chi-square
   if (length(expected[expected<=5])){
      fisher_test_result=fisher.test(x=marker_data[,i], y=target, alternative = "two.sided")
      cont_table_results[i,"cont_table_pvalue"]=fisher_test_result$p.value
      cont_table_results[i,"test"]="Fisher"
   }else{
   #Enter results into array for print out
   cont_table_results[i,"cont_table_pvalue"]=chisquare_result$p.value
   cont_table_results[i,"test"]="Pearson_Chi"
   }
}

#Correct p-values
cont_pvalues=as.numeric(cont_table_results[,"cont_table_pvalue"])
cont_pvalues_adj=mt.rawp2adjp(cont_pvalues, proc=c("Bonferroni","BH"))
cont_pvalues_adj_orig_order=cont_pvalues_adj$adjp[order(cont_pvalues_adj$index),]
cont_table_results=cbind(cont_table_results, cont_pvalues_adj_orig_order[,2:3])

#Write complete results to file
write.table (cont_table_results, sep="\t", file=outfile)

#Unused code
#reduce number of decimal places
#pvalue_results_formatted=format(pvalue_results, digits=5, scientific=FALSE)

