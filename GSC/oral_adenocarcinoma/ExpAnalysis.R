library(multtest)
library(outliers)
library(affy)

print("Loading data")
#Load ENSG info (gene names, status, and type)
setwd("/home/obig/Projects/oral_adenocarcinoma/solexa_data/maq_analysis")
genefile="ENS_v52_ENSG_GeneName_Status_Biotype.txt"
gene_map=read.table(genefile, header = TRUE, na.strings = "NA", sep="\t", row.names=1, as.is=TRUE)

setwd("/home/obig/Projects/oral_adenocarcinoma/solexa_data/maq_analysis/analysis/v3")
datafile="AllLibraries.v3.gene.summary"

outfile="AllLibraries.v3.gene.summary.Ranalysis.txt"
outfile2="AllLibraries.v3.gene.summary.Ranalysis.short.txt"
outfile2B="AllLibraries.v3.gene.summary.Ranalysis.medium.txt"
outfile3="AllLibraries.v3.gene.summary.Ranalysis.interest.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=TRUE)

#Get gene names and other info using ID
genes=data[,1]
gene_lengths=data[,5]
gene_lookup=function(x){
 gene_name=gene_map[x,1]
 gene_status=gene_map[x,2]
 gene_type=gene_map[x,3]
 result=matrix(c(gene_name,gene_status,gene_type), nr=1)
 return(result)
}

print ("Mapping gene IDs to name, status and type")
gene_info=t(apply(as.array(genes),1,gene_lookup))
colnames(gene_info)=c("gene_name","gene_status","gene_type")

print ("Normalizing data")
#Determine smallest 'KNOWN' protein coding gene for gene length normalization
min_gene_length=min(gene_lengths[gene_info[,"gene_status"]=="KNOWN" & gene_info[,"gene_type"]=="protein_coding"])

#Get library data
#lib_data=data[,6:35] #AllLibraries.gene.summary.2
#lib_data=data[,6:38] #AllLibraries.gene.summary.3
lib_data=data[,6:58] #AllLibraries.v3.gene.summary
min_lib_size=min(colSums(lib_data)) #Determine smallest raw library size for normalization purposes

#Put aside raw values for skin lung and blood
skin_met_raw=lib_data[,45]
norm_blood_raw=lib_data[,13]
lung_met_raw=lib_data[,14]

#calculate % of zero-count data for each gene (before adding 1 or normalizing)
perc_zero_fun=function(x){length(x[x<=0])/length(x)}
lib_data_perc_zero=apply(lib_data, 1, perc_zero_fun)

#Add 1 to all data to avoid issues with FC and Fisher p-value calculations
lib_data_plus1=lib_data+1

#Normalize for gene size
gene_length_matrix=matrix(gene_lengths, nrow=length(gene_lengths), ncol=length(lib_data))
lib_data_gene_norm=lib_data_plus1/gene_length_matrix
lib_data_norm=lib_data_gene_norm*min_gene_length

#Normalize for library size: (mapped bases/total mapped bases)*min_lib_size (or arbitrary 1,000,000,000)
lib_data_norm_totals=colSums(lib_data_norm)
lib_data_norm2=lib_data_norm
for (i in 1:length(lib_data_norm_totals)){
 #lib_data_norm2[,i]=(lib_data_norm2[,i]/lib_data_norm_totals[i])*1000000000
 lib_data_norm2[,i]=(lib_data_norm2[,i]/lib_data_norm_totals[i])*min_lib_size
}
#Save row and column names for reapplication after quantiles normalization
cols=colnames(lib_data_norm2)
rows=rownames(lib_data_norm2)

#Quantiles normalize to make distributions more similar
lib_data_norm2=as.data.frame(normalize.quantiles(as.matrix(lib_data_norm2)))
colnames(lib_data_norm2)=cols
rownames(lib_data_norm2)=rows

#Extract normalized compendium and tumor_of_interest data
#compendium=lib_data_norm2[,c(1:21,23:29)]
#skin_met=lib_data_norm2[,30]
compendium=lib_data_norm2[,c(1:12,15:44,46:53)]
skin_met=lib_data_norm2[,45]
norm_blood=lib_data_norm2[,13]
lung_met=lib_data_norm2[,14]

#Issue: what data to use for Fisher p-value calculations?
#1) Use raw values plus 1 (PROBLEM: values to large for Fisher test)
#2) Use gene normalized values plus 1 (PROBLEM: Artificially punishes zero values because Fisher can't use fractions. 
### For gene of size 1000bp if libA=10,000 and LibB=0: 10,000+1/1000=10.001 vs 0+1/1000=0.001.
### Fisher sees this as 10 vs 0 or 1. A much smaller non-significant FC of ~10 instead of real situation of ~10000
#3) Use gene normalized values rescaled so that smallest value is 1 to eliminate all fractions (PROBLEM: values to large for Fisher test)
###Can use if only Chi-square test performed?
##skin_met_norm=lib_data_norm[,45]
#norm_blood_norm=lib_data_norm[,13]
#lung_met_norm=lib_data_norm[,14]
#scaling_factor=1/min(c(skin_met_norm,norm_blood_norm,lung_met_norm))
#skin_met_norm1=skin_met_norm*scaling_factor
#norm_blood_norm1=norm_blood_norm*scaling_factor
#lung_met_norm1=lung_met_norm*scaling_factor

#4) Use gene normalized data without scaling factor
#skin_met_norm1=lib_data_gene_norm[,45]
#norm_blood_norm1=lib_data_gene_norm[,13]
#lung_met_norm1=lib_data_gene_norm[,14]

#4B) Use gene normalized data with scaling factor
#skin_met_norm1=lib_data_norm[,45]+1
#norm_blood_norm1=lib_data_norm[,13]+1
#lung_met_norm1=lib_data_norm[,14]+1

#5) Use gene and library normalized values. Add 1 again to eliminate fractions (PROBLEM: Not sure if library size normalization is proper/fair before Fisher calculation
skin_met_norm1=lib_data_norm2[,45]+1
norm_blood_norm1=lib_data_norm2[,13]+1
lung_met_norm1=lib_data_norm2[,14]+1

#Create compendium plus normal blood and patient met to identify outliers
compendium_lung_norm=lib_data_norm2[,c(14,1:13,15:44,46:53)]
compendium_skin_norm=lib_data_norm2[,c(45,1:13,15:44,46:53)]

print ("Calculating basic summary stats")
#Calculate basic summary stats (means, fold change, stdev, etc)
compendium_means=apply (compendium, 1, mean)
compendium_sds=apply(compendium, 1, sd)
skin_comp_sdds=(skin_met-compendium_means)/compendium_sds
skin_exp_rank=(rank(skin_met, ties.method="min")/length(skin_met))*100
lung_comp_sdds=(lung_met-compendium_means)/compendium_sds
lung_exp_rank=(rank(lung_met, ties.method="min")/length(lung_met))*100

lung_comp_fcs=lung_met/compendium_means
lung_comp_fcs[which(lung_comp_fcs<1 & lung_comp_fcs>0)]=-1/lung_comp_fcs[which(lung_comp_fcs<1 & lung_comp_fcs>0)]
skin_comp_fcs=skin_met/compendium_means
skin_comp_fcs[which(skin_comp_fcs<1 & skin_comp_fcs>0)]=-1/skin_comp_fcs[which(skin_comp_fcs<1 & skin_comp_fcs>0)]
blood_comp_fcs=norm_blood/compendium_means
blood_comp_fcs[which(blood_comp_fcs<1 & blood_comp_fcs>0)]=-1/blood_comp_fcs[which(blood_comp_fcs<1 & blood_comp_fcs>0)]
lung_blood_fcs=lung_met/norm_blood
lung_blood_fcs[which(lung_blood_fcs<1 & lung_blood_fcs>0)]=-1/lung_blood_fcs[which(lung_blood_fcs<1 & lung_blood_fcs>0)]
skin_lung_fcs=skin_met/lung_met
skin_lung_fcs[which(skin_lung_fcs<1 & skin_lung_fcs>0)]=-1/skin_lung_fcs[which(skin_lung_fcs<1 & skin_lung_fcs>0)]
skin_blood_fcs=skin_met/norm_blood
skin_blood_fcs[which(skin_blood_fcs<1 & skin_blood_fcs>0)]=-1/skin_blood_fcs[which(skin_blood_fcs<1 & skin_blood_fcs>0)]

#Create booleans for which data pass filters
#tumor expression greater than 95th percentile
skin_met_top=vector(length=length(skin_met))
skin_met_top[which(skin_met>=quantile(skin_met,0.95))]=TRUE
lung_met_top=vector(length=length(lung_met))
lung_met_top[which(lung_met>=quantile(lung_met,0.95))]=TRUE

#tumor expression less than 50th percentile (less stringent because data is very zero-biased)
skin_met_bottom=vector(length=length(skin_met))
skin_met_bottom[which(skin_met<=quantile(skin_met,0.5))]=TRUE
lung_met_bottom=vector(length=length(lung_met))
lung_met_bottom[which(lung_met<=quantile(lung_met,0.5))]=TRUE

#compendium mean greater than 1
compendium_mean_gt1=vector(length=length(compendium_means))
compendium_mean_gt1[which(compendium_means>=100)]=TRUE

#percent zero values less than 80% or 95%?
perc_zero_filt=vector(length=length(lib_data_perc_zero))
perc_zero_filt[which(lib_data_perc_zero<0.80)]=TRUE

print ("Calculating outlier p-values")
#calculate outlier p-value for each gene
#The patient of interest must be at index 1. We are only interested in p-values for this datapoint
#Define a function to determine outlier p-values
outlier_fun=function(x){
 outlier_scores=scores(x,type="t")
 outlier_pvalues=scores(x, prob=1,type="t") #Uses t-statistic, Z and chisq seem to produce similar results
 patient_outlier_score=outlier_scores[1]
 patient_outlier_pvalue=outlier_pvalues[1] #Return only the p-value for sample of interest
 #The following is necessary because a large positive score is assigned a "pvalue"~1
 #Whereas an equally large -ve score is assigned a "p-value"~0
 if(is.na(patient_outlier_score)==FALSE){if(patient_outlier_score>0){patient_outlier_pvalue=1-patient_outlier_pvalue}}
 return(patient_outlier_pvalue)
}

#Determine outlier p-values for all genes for both skin and lung mets
lung_outlier_pvalues=apply(compendium_lung_norm, 1, outlier_fun)
skin_outlier_pvalues=apply(compendium_skin_norm, 1, outlier_fun)

#Correct pvalues
#skin met
outlier_pvalues_adj=mt.rawp2adjp(skin_outlier_pvalues, proc=c("Bonferroni","BH"))
outlier_pvalues_adj_orig_order_skin=outlier_pvalues_adj$adjp[order(outlier_pvalues_adj$index),]
colnames(outlier_pvalues_adj_orig_order_skin)=c("skin_comp_outlier_rawp","skin_comp_outlier_Bonferroni","skin_comp_outlier_BH")

#lung met
outlier_pvalues_adj2=mt.rawp2adjp(lung_outlier_pvalues, proc=c("Bonferroni","BH"))
outlier_pvalues_adj_orig_order_lung=outlier_pvalues_adj2$adjp[order(outlier_pvalues_adj2$index),]
colnames(outlier_pvalues_adj_orig_order_lung)=c("lung_comp_outlier_rawp","lung_comp_outlier_Bonferroni","lung_comp_outlier_BH")

#Define function to calculate fisher exact statistic
fisher_fun=function(x){
 cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
 fisher_test_result=fisher.test(cont_table, alternative = "two.sided")
 result=matrix(c(fisher_test_result$p.value, fisher_test_result$estimate[[1]]), nr=1)
 return(result)
}

#Define function to calculate fisher or chisq depending on expected counts in contingency table
#Problem with using Chi-Sq instead of Fisher. An approximation with precision limited to p = 1E-16. 
#But, seems to run much faster. Especially for large numbers.
cont_table_fun=function(x){
 cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
 chisquare_result=chisq.test(x=cont_table)
 expected=chisquare_result$expected
 statistic=chisquare_result$statistic[[1]]
 pvalue=chisquare_result$p.value
 if (length(expected[expected<=5])){
  fisher_test_result=fisher.test(x=cont_table, alternative = "two.sided")
  pvalue=fisher_test_result$p.value
  statistic=fisher_test_result$estimate[[1]]
 }
 result=matrix(c(pvalue, statistic), nr=1)
 return(result)
}

#Define a function to calculate chisq p-value only. 
#Given the total counts, if expected counts are low enough to warrant Fisher (i.e., <5), gene is most likely not interesting 
#Assign NA in such cases
chisq_fun=function(x){
 cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
 chisquare_result=chisq.test(x=cont_table)
 expected=chisquare_result$expected
 statistic=chisquare_result$statistic[[1]]
 pvalue=chisquare_result$p.value
 if (length(expected[expected<=5])){
  pvalue=NA
  statistic=NA
 }
 result=matrix(c(pvalue, statistic), nr=1)
 return(result)
}

print("Calculating Fisher p-values")
#Create an array with values for each gene to be analysed by fisher.test
#Use the Normalized coverage values for each library to compare the new skin_met to lung_met, or skin_met to norm_blood
#For the two libraries we will compare the average coverage for geneA versus the sum of average coverage for all other genes
fisher_data_skin_lung=cbind(skin_met_norm1, sum(skin_met_norm1)-skin_met_norm1, lung_met_norm1, sum(lung_met_norm1)-lung_met_norm1)
fisher_data_skin_blood=cbind(skin_met_norm1, sum(skin_met_norm1)-skin_met_norm1, norm_blood_norm1, sum(norm_blood_norm1)-norm_blood_norm1)
fisher_data_lung_blood=cbind(lung_met_norm1, sum(lung_met_norm1)-lung_met_norm1, norm_blood_norm1, sum(norm_blood_norm1)-norm_blood_norm1)

#Apply the fisher function to all rows of the data for the appropriate columns
fisher_results_skin_lung=t(apply(fisher_data_skin_lung, 1, fisher_fun)) #skin_met vs lung_met
fisher_results_skin_blood=t(apply(fisher_data_skin_blood, 1, fisher_fun)) #skin_met vs norm_blood
fisher_results_lung_blood=t(apply(fisher_data_lung_blood, 1, fisher_fun)) #lung_met vs norm_blood

#Correct p-values
fisher_pvalues=as.numeric(fisher_results_skin_lung[,1])
fisher_pvalues_adj=mt.rawp2adjp(fisher_pvalues, proc=c("Bonferroni","BH"))
fisher_pvalues_adj_orig_order=fisher_pvalues_adj$adjp[order(fisher_pvalues_adj$index),]
fisher_results_skin_lung=cbind(fisher_results_skin_lung, fisher_pvalues_adj_orig_order)

fisher_pvalues2=as.numeric(fisher_results_skin_blood[,1])
fisher_pvalues_adj2=mt.rawp2adjp(fisher_pvalues2, proc=c("Bonferroni","BH"))
fisher_pvalues_adj_orig_order2=fisher_pvalues_adj2$adjp[order(fisher_pvalues_adj2$index),]
fisher_results_skin_blood=cbind(fisher_results_skin_blood, fisher_pvalues_adj_orig_order2)

fisher_pvalues3=as.numeric(fisher_results_lung_blood[,1])
fisher_pvalues_adj3=mt.rawp2adjp(fisher_pvalues3, proc=c("Bonferroni","BH"))
fisher_pvalues_adj_orig_order3=fisher_pvalues_adj3$adjp[order(fisher_pvalues_adj3$index),]
fisher_results_lung_blood=cbind(fisher_results_lung_blood, fisher_pvalues_adj_orig_order3)

#Give some column names to the output
colnames(fisher_results_skin_lung)=c("skin_lung_fisher_pvalue", "skin_lung_fisher_odds_ratio","skin_lung_fisher_rawp", "skin_lung_fisher_Bonferroni", "skin_lung_fisher_BH")
colnames(fisher_results_skin_blood)=c("skin_blood_fisher_pvalue", "skin_blood_fisher_odds_ratio","skin_blood_fisher_rawp", "skin_blood_fisher_Bonferroni", "skin_blood_fisher_BH")
colnames(fisher_results_lung_blood)=c("lung_blood_fisher_pvalue", "lung_blood_fisher_odds_ratio","lung_blood_fisher_rawp", "lung_blood_fisher_Bonferroni", "lung_blood_fisher_BH")

print("Summarizing data and writing diff exp gene files")
#Combine all into summary and write to file.
summary=cbind(genes,gene_info,gene_lengths,norm_blood,lung_met,skin_met,compendium_means,lung_comp_fcs,skin_comp_fcs,blood_comp_fcs,skin_lung_fcs,lung_blood_fcs,skin_blood_fcs,lung_comp_sdds,skin_comp_sdds,outlier_pvalues_adj_orig_order_lung,outlier_pvalues_adj_orig_order_skin,fisher_results_lung_blood,fisher_results_skin_blood,fisher_results_skin_lung,lung_exp_rank,skin_exp_rank,lung_met_top,lung_met_bottom,skin_met_top,skin_met_bottom,compendium_mean_gt1,perc_zero_filt,norm_blood_raw,lung_met_raw,skin_met_raw,compendium)
rownames(summary)=genes
write.table (summary, sep="\t", row.names=FALSE, file=outfile)

#Create shorter version of file
summary2=summary[,c("genes","gene_name","compendium_means","norm_blood","lung_met","skin_met","blood_comp_fcs","lung_comp_fcs","skin_comp_fcs","lung_blood_fcs","skin_blood_fcs","skin_lung_fcs","lung_comp_outlier_BH","skin_comp_outlier_BH","skin_lung_fisher_BH","lung_exp_rank","skin_exp_rank","norm_blood_raw","lung_met_raw","skin_met_raw")]
rownames(summary2)=genes
write.table (summary2, sep="\t", row.names=FALSE, file=outfile2)

summary2B=summary[,c("genes","gene_name","gene_status","gene_type","gene_lengths","norm_blood","lung_met","skin_met","compendium_means","lung_comp_fcs","skin_comp_fcs","blood_comp_fcs","skin_lung_fcs","lung_blood_fcs","skin_blood_fcs","lung_comp_sdds","skin_comp_sdds","lung_comp_outlier_rawp","lung_comp_outlier_BH","skin_comp_outlier_rawp","skin_comp_outlier_BH","lung_blood_fisher_BH","skin_blood_fisher_BH","skin_lung_fisher_BH","skin_lung_fisher_Bonferroni","lung_exp_rank","skin_exp_rank","lung_met_top","lung_met_bottom","skin_met_top","skin_met_bottom","compendium_mean_gt1","perc_zero_filt")]
rownames(summary2B)=genes
write.table (summary2B, sep="\t", row.names=FALSE, file=outfile2B)

#Create output file with just genes of interest
#summary3=summary2[c("ENSG00000187840","ENSG00000142208","ENSG00000105221","ENSG00000002330","ENSG00000157764","ENSG00000157613","ENSG00000182578","ENSG00000138798","ENSG00000058085","ENSG00000100030","ENSG00000102882","ENSG00000169032","ENSG00000167034","ENSG00000149269","ENSG00000180370","ENSG00000137843","ENSG00000197461","ENSG00000100311","ENSG00000073756","ENSG00000173801","ENSG00000197943","ENSG00000075651","ENSG00000112033","ENSG00000171862","ENSG00000128340","ENSG00000169750","ENSG00000132155","ENSG00000006451","ENSG00000116473","ENSG00000127314","ENSG00000174775","ENSG00000213281","ENSG00000139687","ENSG00000165731","ENSG00000146648","ENSG00000198793","ENSG00000151247","ENSG00000135930","ENSG00000108443"),]

genelist1=c("ENSG00000142208","ENSG00000105221","ENSG00000117020","ENSG00000140992","ENSG00000171862","ENSG00000121879","ENSG00000051382","ENSG0000071608","ENSG00000105851","ENSG00000198793","ENSG00000164327","ENSG00000141564","ENSG00000063046","ENSG00000151247","ENSG00000175766","ENSG00000135930","ENSG00000163412","ENSG00000187840","ENSG00000148730","ENSG00000100784","ENSG00000108443","ENSG00000175634","ENSG00000117676","ENSG00000071242","ENSG00000177189","ENSG00000162302","ENSG00000072133","ENSG00000136643","ENSG00000138798","ENSG00000146648","ENSG00000197461","ENSG00000100311","ENSG00000145431","ENSG00000170962","ENSG00000134853","ENSG00000113721","ENSG00000165731","ENSG00000109458","ENSG00000033327")
genelist2=c("ENSG00000177885","ENSG00000160691","ENSG00000129946","ENSG00000148082","ENSG00000185634","ENSG00000115904","ENSG00000100485","ENSG00000174775","ENSG00000133703","ENSG00000213281","ENSG00000157764","ENSG00000132155","ENSG00000169032","ENSG00000126934","ENSG00000100030","ENSG00000102882","ENSG00000100784","ENSG00000108443","ENSG00000175634","ENSG00000167034","ENSG00000141510","ENSG00000139687","ENSG00000080031","ENSG00000101773")
genelist=c(genelist1,genelist2)
summary3=summary2[genelist,]
write.table (summary3, sep="\t", row.names=FALSE,file=outfile3)


#Create lung/skin up/down-regulated (vs compendium) and skin_vs_lung gene lists
#First, apply filters: (compendium_mean>1 or zero_perc filter) and percentile filters
#filtered_comp=summary2B[summary2B[,"compendium_mean_gt1"]==TRUE,]
filtered_comp=summary2B[summary2B[,"perc_zero_filt"]==TRUE,]

#Filter out UNKNOWN and non-protein_coding genes
filtered_comp=filtered_comp[filtered_comp[,"gene_status"]=="KNOWN" & filtered_comp[,"gene_type"]=="protein_coding",]

#Determine outlier p-value thresholds to use for down-regulation that are 'equivalent' to up-regulation threshold of 0.05
#Use raw p-values because these have more precision, less ties than BH values
#Lung
#filtered_lung=filtered_comp
#pvalues=as.numeric(filtered_lung[,"lung_comp_outlier_rawp"])
#pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
#pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
#filtered_lung[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_lung_UP=filtered_lung[filtered_lung[,"lung_comp_fcs"]>0,]
#filtered_lung_DOWN=filtered_lung[filtered_lung[,"lung_comp_fcs"]<0,]
#filtered_lung_sig_UP=filtered_lung_UP[filtered_lung_UP[,"lung_comp_outlier_BH"]<0.05,]
#lung_UP_count=length(filtered_lung_sig_UP[,1])
#down_outlier_threshold_lung=sort(filtered_lung_DOWN[,"lung_comp_outlier_rawp"])[lung_UP_count]
down_outlier_threshold_lung=0.10

#Skin
#filtered_skin=filtered_comp
#pvalues=as.numeric(filtered_skin[,"skin_comp_outlier_rawp"])
#pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
#pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
#filtered_skin[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_skin_UP=filtered_skin[filtered_skin[,"skin_comp_fcs"]>0,]
#filtered_skin_DOWN=filtered_skin[filtered_skin[,"skin_comp_fcs"]<0,]
#filtered_skin_sig_UP=filtered_skin_UP[filtered_skin_UP[,"skin_comp_outlier_BH"]<0.05,]
#skin_UP_count=length(filtered_skin_sig_UP[,1])
#down_outlier_threshold_skin=sort(filtered_skin_DOWN[,"skin_comp_outlier_rawp"])[skin_UP_count]
down_outlier_threshold_skin=0.10

#Lung - up
#with expression percentile filter
filtered_lung_top=filtered_comp[filtered_comp[,"lung_met_top"]==TRUE,]
pvalues=as.numeric(filtered_lung_top[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung_top[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_lung_sig=filtered_lung_top[filtered_lung_top[,"lung_comp_outlier_BH"]<0.05,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]>2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]>1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_up_top.txt")

#without expression percentile filter
filtered_lung=filtered_comp
pvalues=as.numeric(filtered_lung[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_lung_sig=filtered_lung[filtered_lung[,"lung_comp_outlier_BH"]<0.05,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]>2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]>1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_up.txt")

#without expression percentile filter, and with FC filter (before MTC)
filtered_lung=filtered_comp
filtered_lung_FC=filtered_lung[filtered_lung[,"lung_comp_fcs"]>0,]
pvalues=as.numeric(filtered_lung_FC[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung_FC[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_lung_sig=filtered_lung_FC[filtered_lung_FC[,"lung_comp_outlier_BH"]<0.05,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]>2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]>1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_up_FCprefilt.txt")

#Lung - down
#with expression percentile filter
filtered_lung_bottom=filtered_comp[filtered_comp[,"lung_met_bottom"]==TRUE,]
pvalues=as.numeric(filtered_lung_bottom[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung_bottom[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_lung_sig=filtered_lung_bottom[filtered_lung_bottom[,"lung_comp_outlier_BH"]<0.10,]
filtered_lung_sig=filtered_lung_bottom[filtered_lung_bottom[,"lung_comp_outlier_rawp"]<down_outlier_threshold_lung,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]< -2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]< -1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_down_bottom.txt")

#without expression percentile filter
filtered_lung=filtered_comp
pvalues=as.numeric(filtered_lung[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_lung_sig=filtered_lung[filtered_lung[,"lung_comp_outlier_BH"]<0.10,]
filtered_lung_sig=filtered_lung[filtered_lung[,"lung_comp_outlier_rawp"]<down_outlier_threshold_lung,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]< -2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]< -1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_down.txt")

#without expression percentile filter, and with FC filter (before MTC)
filtered_lung=filtered_comp
filtered_lung_FC=filtered_lung[filtered_lung[,"lung_comp_fcs"]<0,]
pvalues=as.numeric(filtered_lung_FC[,"lung_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_lung_FC[,"lung_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_lung_sig=filtered_lung_FC[filtered_lung_FC[,"lung_comp_outlier_BH"]<0.10,]
filtered_lung_sig=filtered_lung_FC[filtered_lung_FC[,"lung_comp_outlier_rawp"]<down_outlier_threshold_lung,]
filtered_lung_sig=filtered_lung_sig[!is.na(filtered_lung_sig[,"lung_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_lung_sig2=filtered_lung_sig[filtered_lung_sig[,"lung_blood_fisher_BH"]<0.05,]
filtered_lung_sig_FC=filtered_lung_sig2[filtered_lung_sig2[,"lung_comp_fcs"]< -2,]
filtered_lung_sig_FC2=filtered_lung_sig_FC[filtered_lung_sig_FC[,"lung_blood_fcs"]< -1.5,]
write.table (filtered_lung_sig_FC2, sep="\t", row.names=FALSE, file="lung_down_FCprefilt.txt")

#Skin - up
#with expression percentile filter
filtered_skin_top=filtered_comp[filtered_comp[,"skin_met_top"]==TRUE,]
pvalues=as.numeric(filtered_skin_top[,"skin_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_skin_top[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_skin_sig=filtered_skin_top[filtered_skin_top[,"skin_comp_outlier_BH"]<0.05,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]>2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]>1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_up_top.txt")

#without expression percentile filter
filtered_skin=filtered_comp
pvalues=as.numeric(filtered_skin[,"skin_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_skin[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_skin_sig=filtered_skin[filtered_skin[,"skin_comp_outlier_BH"]<0.05,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]>2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]>1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_up.txt")

##without expression percentile filter, and with FC filter (before MTC)
filtered_skin=filtered_comp
filtered_skin_FC=filtered_skin[filtered_skin[,"skin_comp_fcs"]>0,]
pvalues=as.numeric(filtered_skin_FC[,"skin_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_skin_FC[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
filtered_skin_sig=filtered_skin_FC[filtered_skin_FC[,"skin_comp_outlier_BH"]<0.05,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]>2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]>1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_up_FCprefilt.txt")

#Skin - down
#with expression percentile filter
filtered_skin_bottom=filtered_comp[filtered_comp[,"skin_met_bottom"]==TRUE,]
#pvalues=as.numeric(filtered_skin_bottom[,"skin_comp_outlier_rawp"])
#pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
#pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
#filtered_skin_bottom[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_skin_sig=filtered_skin_bottom[filtered_skin_bottom[,"skin_comp_outlier_BH"]<0.10,]
filtered_skin_sig=filtered_skin_bottom[filtered_skin_bottom[,"skin_comp_outlier_rawp"]<down_outlier_threshold_skin,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]< -2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]< -1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_down_bottom.txt")

#without expression percentile filter
filtered_skin=filtered_comp
pvalues=as.numeric(filtered_skin[,"skin_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_skin[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_skin_sig=filtered_skin[filtered_skin[,"skin_comp_outlier_BH"]<0.10,]
filtered_skin_sig=filtered_skin[filtered_skin[,"skin_comp_outlier_rawp"]<down_outlier_threshold_skin,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]< -2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]< -1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_down.txt")

#without expression percentile filter, and with FC filter (before MTC)
filtered_skin=filtered_comp
filtered_skin_FC=filtered_skin[filtered_skin[,"skin_comp_fcs"]<0,]
pvalues=as.numeric(filtered_skin_FC[,"skin_comp_outlier_rawp"])
pvalues_adj=mt.rawp2adjp(pvalues, proc="BH")
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
filtered_skin_FC[,"skin_comp_outlier_BH"]=pvalues_adj_orig_order[,"BH"]
#filtered_skin_sig=filtered_skin_FC[filtered_skin_FC[,"skin_comp_outlier_BH"]<0.10,]
filtered_skin_sig=filtered_skin_FC[filtered_skin_FC[,"skin_comp_outlier_rawp"]<down_outlier_threshold_skin,]
filtered_skin_sig=filtered_skin_sig[!is.na(filtered_skin_sig[,"skin_blood_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_sig2=filtered_skin_sig[filtered_skin_sig[,"skin_blood_fisher_BH"]<0.05,]
filtered_skin_sig_FC=filtered_skin_sig2[filtered_skin_sig2[,"skin_comp_fcs"]< -2,]
filtered_skin_sig_FC2=filtered_skin_sig_FC[filtered_skin_sig_FC[,"skin_blood_fcs"]< -1.5,]
write.table (filtered_skin_sig_FC2, sep="\t", row.names=FALSE, file="skin_down_FCprefilt.txt")

#Skin vs Lung - up/down - with expression filter
filtered_skin_lung_top=filtered_comp[filtered_comp[,"skin_met_top"]==TRUE | filtered_comp[,"lung_met_top"]==TRUE,]
filtered_skin_lung_top=filtered_skin_lung_top[!is.na(filtered_skin_lung_top[,"skin_lung_fisher_BH"]),] #Filter out Fisher NAs
filtered_skin_lung_sig=filtered_skin_lung_top[filtered_skin_lung_top[,"skin_lung_fisher_BH"]<0.05,]
filtered_skin_lung_sig_FC_UP=filtered_skin_lung_sig[filtered_skin_lung_sig[,"skin_lung_fcs"]>2,]
filtered_skin_lung_sig_FC_DOWN=filtered_skin_lung_sig[filtered_skin_lung_sig[,"skin_lung_fcs"]< -2,]
write.table (filtered_skin_lung_sig_FC_UP, sep="\t", row.names=FALSE, file="skin_vs_lung_up_top.txt")
write.table (filtered_skin_lung_sig_FC_DOWN, sep="\t", row.names=FALSE, file="skin_vs_lung_down_top.txt")

#Skin vs Lung - up/down - without expression filter
filtered_skin_lung=filtered_comp
#filtered_skin_lung=filtered_skin_lung[!is.na(filtered_skin_lung[,"skin_lung_fisher_BH"]),] #Filter out Fisher NAs
#filtered_skin_lung_sig=filtered_skin_lung[filtered_skin_lung[,"skin_lung_fisher_BH"]<0.05,] #use more stringent p-value?
filtered_skin_lung=filtered_skin_lung[!is.na(filtered_skin_lung[,"skin_lung_fisher_Bonferroni"]),] #Filter out Fisher NAs
filtered_skin_lung_sig=filtered_skin_lung[filtered_skin_lung[,"skin_lung_fisher_Bonferroni"]<0.001,] #use more stringent p-value?
filtered_skin_lung_sig_FC_UP=filtered_skin_lung_sig[filtered_skin_lung_sig[,"skin_lung_fcs"]>2,]
filtered_skin_lung_sig_FC_DOWN=filtered_skin_lung_sig[filtered_skin_lung_sig[,"skin_lung_fcs"]< -2,]
write.table (filtered_skin_lung_sig_FC_UP, sep="\t", row.names=FALSE, file="skin_vs_lung_up.txt")
write.table (filtered_skin_lung_sig_FC_DOWN, sep="\t", row.names=FALSE, file="skin_vs_lung_down.txt")

















#Unused code

#summary2=cbind(gene_lengths,norm_blood,lung_met,skin_met,compendium_means,lung_comp_fcs,skin_comp_fcs,blood_comp_fcs,skin_lung_fcs,lung_blood_fcs,skin_blood_fcs,lung_comp_sdds,skin_comp_sdds,wilcox_pvalues_adj_orig_order2[,"BH"],wilcox_pvalues_adj_orig_order[,"BH"],outlier_pvalues_adj_orig_order_lung[,c("rawp","BH")],outlier_pvalues_adj_orig_order_skin[,c("rawp","BH")],fisher_results_lung_blood[,"BH"],fisher_results_skin_blood[,"BH"],fisher_results_skin_lung[,"BH"],lung_exp_rank,skin_exp_rank,lung_met_top5,lung_met_bottom5,skin_met_top5,skin_met_bottom5,compendium_mean_gt1)
#rownames(summary2)=genes
#colnames(summary2)=c("gene_length","norm_blood","lung_met","skin_met","compendium_mean","lung_comp_FC","skin_comp_FC","blood_comp_FC","skin_lung_FC","lung_blood_FC","skin_blood_FC","lung_comp_SDD","skin_comp_SDD","lung_comp_wilcox_BH","skin_comp_wilcox_BH","lung_comp_outlier_rawp","lung_comp_outlier_BH","skin_comp_outlier_rawp","skin_comp_outlier_BH","lung_blood_fisher_BH","skin_blood_fisher_BH","skin_lung_fisher_BH","lung_exp_rank","skin_exp_rank","lung_top5","lung_bottom5","skin_top5","skin_bottom5","compendium_mean_gt1")
#write.table (summary2, sep="\t", file=outfile2)

#calculate outlier p-value for each gene
#Outliers are determined one at a time (e.g.,value most different from mean)
#Remove or replace (e.g., with mean or median) each outlier and determine the 'next' outlier
#The problem with replacing is that sometimes the replaced value later becomes an 'outlier' again.
#The problem with removing is that sometimes the second most extreme outlier has a more significant p-value than the most extreme
#The patient of interest is at index 1. We are only interested in p-values for this datapoint

#Create a datastructure for outlier p-values
#libraries=colnames(compendium_lung_norm)
#result=array(0,dimnames=list(libraries,c("pvalue", "value")), dim=c(length(libraries),2))

#values=as.numeric(compendium_lung_norm[1,])
#x=values
#num_tests=length(x)

#for (i in 1:num_tests){
# outlier=outlier(x)
# outlier_index=which(values==outlier)
# outlier_test=grubbs.test(x)
# result[outlier_index,"value"]=outlier
# result[outlier_index,"pvalue"]=outlier_test$p.value
# x=rm.outlier(x,fill=FALSE) #Replace outlier with mean or median (which to use probably depends on what outlier test makes use of)
#}


#Perform a two-tailed wilcox test where mu=tumor_of_interest.
##skin versus compendium
#wilcox_u_fun = function(x){
##  wilcox.test(x[c(1:21,23:29)], mu=x[30])$p.value
##  wilcox.test(x[c(1:21,26:32)], mu=x[33])$p.value
#  wilcox.test(x[c(1:13,15:44,46:53)], mu=x[45])$p.value
#}
#wilcox_test_pvalues=apply(lib_data_norm2, 1, wilcox_u_fun)

#Correct pvalues
#wilcox_pvalues_adj=mt.rawp2adjp(wilcox_test_pvalues, proc=c("Bonferroni","BH"))
#wilcox_pvalues_adj_orig_order_skin=wilcox_pvalues_adj$adjp[order(wilcox_pvalues_adj$index),]
#colnames(wilcox_pvalues_adj_orig_order_skin)=c("skin_comp_wilcox_rawp","skin_comp_wilcox_Bonferroni","skin_comp_wilcox_BH")

##lung versus compendium
#wilcox_u_fun2 = function(x){
##  wilcox.test(x[c(1:21,26:32)], mu=x[24])$p.value
#  wilcox.test(x[c(1:13,15:44,46:53)], mu=x[14])$p.value
#}
#wilcox_test_pvalues2=apply(lib_data_norm2, 1, wilcox_u_fun2)

#Correct pvalues
#wilcox_pvalues_adj2=mt.rawp2adjp(wilcox_test_pvalues2, proc=c("Bonferroni","BH"))
#wilcox_pvalues_adj_orig_order_lung=wilcox_pvalues_adj2$adjp[order(wilcox_pvalues_adj2$index),]
#colnames(wilcox_pvalues_adj_orig_order_lung)=c("lung_comp_wilcox_rawp","lung_comp_wilcox_Bonferroni","lung_comp_wilcox_BH")


