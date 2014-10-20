library(multtest)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"

#outfile="BenignMalignant_58markers_23JAN08_group1_progtests.2.txt"
outfile="BenignMalignant_58markers_23JAN08_group2_progtests.2.txt"
#outfile="BenignMalignant_58markers_23JAN08_group1_progtests.txt"
#outfile="BenignMalignant_58markers_23JAN08_group2_progtests.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
path=data[,13]
data=data[path<4,] #Exclude M malignancies
path=data[,13]
data=data[path>0,] #Exclude Benign

#Get marker data
#marker_data=data[,36:92] #ungrouped values
#marker_data=data[,93:149] #group 1 values
marker_data=data[,150:206] #group 2 values
#marker_data=data[,c("Galectin_g1","CK19_g1","VEGF_g1","AuroraA_g1","P16_g1","AR_g1","HBME_g1","BCL2_g1","CYCLIND1_g1","Caveolin1_g1","CYCLINE_g1","ECAD_g1","CR3_g1","Clusterin_g1","IGFBP5_g1","P21_g1","IGFBP2_g1","BetaCatenin_g1","HER4_g1","TG_g1","KI67_g1","Caveolin_g1","AuroraC_g1","S100_g1","MRAS_g1","CKIT_g1","HER3_g1","RET_g1","AMFR_g1","MLH1_g1","AAT_g1","TTF1_g1","PGI_g1","HSP27_g1","Syntrophin_g1")] #group 1 values for sig markers from diag tests
#marker_data=data[,c("Galectin_g2","CK19_g2","VEGF_g2","AuroraA_g2","P16_g2","AR_g2","HBME_g2","BCL2_g2","CYCLIND1_g2","Caveolin1_g2","CYCLINE_g2","ECAD_g2","CR3_g2","Clusterin_g2","IGFBP5_g2","P21_g2","IGFBP2_g2","BetaCatenin_g2","HER4_g2","TG_g2","KI67_g2","Caveolin_g2","AuroraC_g2","S100_g2","MRAS_g2","CKIT_g2","HER3_g2","RET_g2","AMFR_g2","MLH1_g2","AAT_g2","TTF1_g2","PGI_g2","HSP27_g2","Syntrophin_g2")] #group 2 values for sig markers from diag tests
marker_names=colnames(marker_data)

#Get categorical prognostic variables
#prognostic_data=data[,c("Sex","age_group","Fam_Hx","Path","LN_Met","ETE","Size_grp","Vasc","Multifoc","Distant_Mets","T","N","M","AMES","RAI","EBRT","T4","Dz_Status")]
prognostic_data=data[,c("Sex","Path","ETE","Vasc","Complete_Resect","T","N","M","AMES")] #Remove duplicate/uninteresting variables to reduce number of tests, only include age and size as continuous.
prognostic_names=colnames(prognostic_data)

#Get continuous prognostic variables
cont_prognostic_data=data[,c("Age","Size_Ca","MACIS")]
cont_prognostic_names=colnames(cont_prognostic_data)

#Make sure categorical data is correctly recognized as categorical
for(i in 1:length(prognostic_data)){
prognostic_data[,i]=as.factor(prognostic_data[,i])
}
for(i in 1:length(marker_data)){
marker_data[,i]=as.factor(marker_data[,i])
}

#Get list of categorical variable tests to be run for results array
marker_vs_prog_var=vector(length=length(marker_names)*length(prognostic_names))
n=0
for (i in 1:length(marker_names)){
 for (j in 1:length(prognostic_names)){
  n=n+1
  marker_vs_prog_var[n]=paste(marker_names[i],"_x_",prognostic_names[j],sep="")
 }
}

#Get list of continuous variable tests to be run for results array
marker_vs_cont_prog_var=vector(length=length(marker_names)*length(cont_prognostic_names))
n=0
for (i in 1:length(marker_names)){
 for (j in 1:length(cont_prognostic_names)){
  n=n+1
  marker_vs_cont_prog_var[n]=paste(marker_names[i],"_x_",cont_prognostic_names[j],sep="")
 }
}

#Combine all test names
marker_vs_all_prog_var=c(marker_vs_prog_var, marker_vs_cont_prog_var)

#Create array to store results for prognostic contingency table stats
table_results = array(0, dimnames = list(marker_vs_all_prog_var, c("pvalue", "test")), dim=c(length(marker_vs_all_prog_var),2))

#Loop through all markers for categorical tests
n=0
for (i in 1:length(marker_names)){
 #loop through all categorical prognostic variables and perform contingency table statistics
 for (j in 1:length(prognostic_names)){
   n=n+1
   #Skip any markers which do not have at least two factors (i.e. if constant)
   if (nlevels(marker_data[,i])<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   #Skip any prognostic variables which do not have at least two factors (i.e. if constant)
   if (nlevels(prognostic_data[,j])<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   
   #Skip any marker/prognostic variable combinations which do not have at least two factors when cases are excluded (i.e. due to missing values)
   ###These can probably replace the two checks done above
   if (nlevels(as.factor(as.vector(marker_data[,i][!is.na(prognostic_data[,j])])))<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   if (nlevels(as.factor(as.vector(prognostic_data[,j][!is.na(marker_data[,i])])))<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }

   #Calculate Pearson Chi-square for all markers vs diagnosis
   chisquare_result=chisq.test(x=marker_data[,i], y=prognostic_data[,j])
   observed=chisquare_result$observed
   expected=chisquare_result$expected
   #If one or more cells in the table have expected count less than or equal to 5 use a Fisher test instead of Pearson Chi-square
   if (length(expected[expected<=5])>0){
      fisher_test_result=fisher.test(x=marker_data[,i], y=prognostic_data[,j], alternative = "two.sided")
      table_results[n,"pvalue"]=fisher_test_result$p.value
      table_results[n,"test"]="Fisher"
   }else{
   #Enter results into array for print out
   table_results[n,"pvalue"]=chisquare_result$p.value
   table_results[n,"test"]="Pearson_Chi"
  }
 }
}

#Loop through all markers for continuous tests
for (i in 1:length(marker_names)){
 #loop through all continuous prognostic variables and perform Mann-Whitney statistics
 for (k in 1:length(cont_prognostic_names)){
   n=n+1
   #Skip any markers which do not have at least two factors (i.e. if constant)
   if (nlevels(marker_data[,i])<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   #Get continuous variable values for two marker groups ("low" and "high")
   low_marker_values=cont_prognostic_data[marker_data[,i]==0,k]
   high_marker_values=cont_prognostic_data[marker_data[,i]==1,k]
   #Make sure that there are some non-NA values to run the statistic on, otherwise it fails
   if (length(low_marker_values)-sum(complete.cases(low_marker_values)==FALSE)<1){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   if (length(high_marker_values)-sum(complete.cases(high_marker_values)==FALSE)<1){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
   }
   #Calculate Mann-Whitney U-test for two marker score groups versus continuous prognostic variables
   wilcox_result=wilcox.test(x=low_marker_values, y=high_marker_values, alternative="two.sided", paired=FALSE)
   table_results[n,"pvalue"]=wilcox_result$p.value
   table_results[n,"test"]="Mann-Whitney"
 }
}

#Remove rows where no p-value could be calculated.  This screws up multtest even though it shouldn't matter
table_results=table_results[!is.na(table_results[,1]),]

#Correct p-values
pvalues=as.numeric(table_results[,"pvalue"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
table_results=cbind(table_results, pvalues_adj_orig_order[,2:3])

#Write complete results to file
write.table (table_results, sep="\t", file=outfile)

#Unused code
#reduce number of decimal places
#pvalue_results_formatted=format(pvalue_results, digits=5, scientific=FALSE)

