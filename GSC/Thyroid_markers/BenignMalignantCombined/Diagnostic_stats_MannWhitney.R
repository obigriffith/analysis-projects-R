library(multtest)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"
outfile="BenignMalignant_58markers_23JAN08_ungrouped_MannWhitney.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
path=data[,13]

#Keep only rows where pathology is 0, 1, 2, or 3. Note we are only excluding M cases not HCC as previously done
data=data[path<4,]

#Get marker data
marker_data=data[,36:92] #ungrouped values

marker_names=colnames(marker_data)

#Get diagnostic target variable and specify as factor/categorical
target=(data[,14])
target[target==0]="Benign"
target[target==1]="Malignant"
target=as.factor(target)

#Create array to store results for diagnostic contingency table stats
wilcox_results = array(0, dimnames = list(marker_names, c("wilcox_pvalue", "W", "Benign_mean_rank", "Malignant_mean_rank", "Direction")), dim=c(length(marker_names),5))

#For each marker perform the Mann-Whitney U test for marker score vs. pathology (Benign vs. Malignant)
for (i in 1:length(marker_data)){
   #Summarize data for marker vs diagnosis
   #Calculate mean ranks for Benign vs Malignant so that direction of change can be determined
   ranked_marker_data=rank(marker_data[,i], na.last="keep")
   benign_mean_rank=mean(ranked_marker_data[target=="Benign"], na.rm=TRUE)
   malig_mean_rank=mean(ranked_marker_data[target=="Malignant"], na.rm=TRUE)
   wilcox_results[i,"Benign_mean_rank"]=benign_mean_rank
   wilcox_results[i,"Malignant_mean_rank"]=malig_mean_rank
   #Determine if marker score has increased or decreased
   if(malig_mean_rank-benign_mean_rank>0){
      direction="UP"
   }
   if(malig_mean_rank-benign_mean_rank<0){
      direction="DOWN"
   }
   if(malig_mean_rank-benign_mean_rank==0){
      direction=NA
   }
   wilcox_results[i,"Direction"]=direction

   #Perform Mann-Whitney U test
   Benign_scores=marker_data[target=="Benign",i]
   Malig_scores=marker_data[target=="Malignant",i]
   wilcox_result=wilcox.test(x=Benign_scores, y=Malig_scores, alternative="two.sided", paired=FALSE)
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

