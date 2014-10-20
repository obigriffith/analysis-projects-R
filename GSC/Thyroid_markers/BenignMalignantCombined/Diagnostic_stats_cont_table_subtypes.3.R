library(multtest)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis/data_files")
datafile="BenignMalignant_58markers_23JAN08_subtype_group_indeterm_RFmarg.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Change dirs for output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis/diag_tests")

outfile="BenignMalignant_58markers_23JAN08_group1_subtypes3_diagtests_fix.txt"
#outfile="BenignMalignant_58markers_23JAN08_group2_subtypes3_diagtests_fix.txt"

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
path=data[,13]

#Keep only rows where pathology is 0, 1, 2, or 3. Note we are only excluding M cases not HCC as previously done
data=data[path<4,]

#Get marker data (Check to make sure columns are correct)
#marker_data=data[,42:98] #ungrouped values
marker_data=data[,99:155] #group 1 values
#marker_data=data[,156:212] #group 2 values

marker_names=colnames(marker_data)

#Get diagnostic target variable and specify as factor/categorical
target=(data[,18])
target=as.factor(target)

#Get rid of data for which no valid subtype exists
marker_data=marker_data[!is.na(target),]
target=target[!is.na(target)]

#Make sure predictor data is correctly recognized as categorical
for(i in 1:length(marker_data)){
marker_data[,i]=as.factor(marker_data[,i])
}

#Create array to store results for diagnostic contingency table stats
cont_table_results = array(0, dimnames = list(marker_names, c("Ben_Neg", "Ben_Pos", "Mal_Neg", "Mal_Pos", "Ben_Pos_perc", "Mal_Pos_perc", "cont_table_pvalue", "test")), dim=c(length(marker_names),8))

for (i in 1:length(marker_data)){
   #Summarize data for marker vs diagnosis
   Benign_scores=marker_data[target=="Foll_Ben",i]
   Malig_scores=marker_data[target=="Foll_Mal",i]
   Benign_scores_noNA=Benign_scores[!is.na(Benign_scores)]
   Malig_scores_noNA=Malig_scores[!is.na(Malig_scores)]
   Ben_Neg=length(Benign_scores_noNA[Benign_scores_noNA==0])
   Ben_Pos=length(Benign_scores_noNA[Benign_scores_noNA==1])
   Mal_Neg=length(Malig_scores_noNA[Malig_scores_noNA==0])
   Mal_Pos=length(Malig_scores_noNA[Malig_scores_noNA==1])
   Ben_Pos_perc=(Ben_Pos/(Ben_Pos+Ben_Neg))*100
   Mal_Pos_perc=(Mal_Pos/(Mal_Pos+Mal_Neg))*100
   cont_table_results[i,"Ben_Neg"]=Ben_Neg
   cont_table_results[i,"Ben_Pos"]=Ben_Pos
   cont_table_results[i,"Mal_Neg"]=Mal_Neg
   cont_table_results[i,"Mal_Pos"]=Mal_Pos
   cont_table_results[i,"Ben_Pos_perc"]=Ben_Pos_perc
   cont_table_results[i,"Mal_Pos_perc"]=Mal_Pos_perc

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

