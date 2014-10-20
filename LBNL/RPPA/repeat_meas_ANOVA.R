library(ez)
library(multtest)

#Use the following to ensure your sum of squares calculations match what you get, for example, in SPSS, or in many common textbooks on ANOVA
options(contrasts=c("contr.sum","contr.poly"))

setwd("C:/Users/Obi/Documents/Projects/HER2_Breast_Cancer/RPPA/R_analysis")

#Load data
datafile="JoeGray480_MedianCentered_179Ab.clean.data.txt"
rawdata=read.table(datafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t")

#Load protein details
protein_metafile="JoeGray480_MedianCentered_179Ab.clean.protein_meta.txt"
protein_meta=read.table(protein_metafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t", as.is=TRUE)
proteins=colnames(protein_meta)
protein_names=protein_meta[1,]

#Load cell line details
lines_metafile="JoeGray480_MedianCentered_179Ab.clean.lines_meta.txt"
lines_meta_raw=read.table(lines_metafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t")

##choose either (1) or (2) and set the results file
#1. All lines
data=rawdata
lines_meta=lines_meta_raw
ANOVA_outfile="synergy_vs_nonsynergy_all_ANOVA.txt"

#2. HER2+ lines only (Amplified only)
data=rawdata[rawdata[,"HER2_status"]=="Amplified",]
lines_meta=lines_meta_raw[lines_meta_raw[,"HER2_status"]=="Amplified",]
ANOVA_outfile="synergy_vs_nonsynergy_HER2amp_ANOVA.txt"

#Perform ANOVA for all proteins
ANOVA_results = array(0, dimnames = list(c("treatment","synergy","treatment*synergy"), proteins), dim=c(3,length(proteins)))

#factors and subject ids are same for each protein
wf1 = as.character(data[,"time_hrs"])
bf1 = as.character(data[,"treatment"])
bf2 = as.character(data[,"Synergy_status"])
sid = as.character(data[,"line_treat_id"])

for (i in 1:length(proteins)){
 #Create dataframe compatible with ezANOVA
 x1 = as.numeric(data[,proteins[i]]) #Abundance data for protein X
 xdf<-data.frame(x1,wf1, bf1, bf2, sid)
 names(xdf)=c("abundance", "time_hrs", "treatment", "synergy", "line_treat_id")

 #Calculate ANOVA for mixed within-between factors
 ezout = ezANOVA(data=xdf, dv=.(abundance), sid=.(line_treat_id), within=.(time_hrs), between=.(treatment, synergy)) 

 #ANOVA Results: ezout$ANOVA
 ANOVA_results["treatment",proteins[i]]=ezout$ANOVA[ezout$ANOVA[,"Effect"]=="treatment","p"]
 ANOVA_results["synergy",proteins[i]]=ezout$ANOVA[ezout$ANOVA[,"Effect"]=="synergy","p"]
 ANOVA_results["treatment*synergy",proteins[i]]=ezout$ANOVA[ezout$ANOVA[,"Effect"]=="treatment:synergy","p"]
}


#CORRECT FOR MULTIPLE TESTING
#Correct p-values
treatment_pvalues=as.numeric(ANOVA_results["treatment",])
treatment_pvalues_adj=mt.rawp2adjp(treatment_pvalues, proc=c("BH"))
treatment_pvalues_adj_orig_order=treatment_pvalues_adj$adjp[order(treatment_pvalues_adj$index),]
ANOVA_results=rbind(ANOVA_results, treatment_pvalues_adj_orig_order[,2])

synergy_pvalues=as.numeric(ANOVA_results["synergy",])
synergy_pvalues_adj=mt.rawp2adjp(synergy_pvalues, proc=c("BH"))
synergy_pvalues_adj_orig_order=synergy_pvalues_adj$adjp[order(synergy_pvalues_adj$index),]
ANOVA_results=rbind(ANOVA_results, synergy_pvalues_adj_orig_order[,2])

treatment_synergy_pvalues=as.numeric(ANOVA_results["treatment*synergy",])
treatment_synergy_pvalues_adj=mt.rawp2adjp(treatment_synergy_pvalues, proc=c("BH"))
treatment_synergy_pvalues_adj_orig_order=treatment_synergy_pvalues_adj$adjp[order(treatment_synergy_pvalues_adj$index),]
ANOVA_results=rbind(ANOVA_results, treatment_synergy_pvalues_adj_orig_order[,2])

#Write results to file
ANOVA_out=t(rbind(protein_names,ANOVA_results))
write.table (ANOVA_out, sep="\t", file=ANOVA_outfile)





#Other data of potential interest
#Mauchly's Test for Sphericity: ezout$Mauchly
#Sphericity Corrections: ezout$Spher

#Calculate Mean, SD, etc for "between" comparisons
ezstats_between=ezStats(data=xdf, dv=.(abundance), between=.(treatment, synergy),  sid = .(line_treat_id))

#Calculate Mean, SD, etc for "within" comparisons
ezstats_within=ezStats(data=xdf, dv=.(abundance), within = .(time_hrs),  sid = .(line_treat_id))

###PLOTTING###
#Plot differences for "between" comparisons (i.e., treatment, synergy_status)
#treatment
#ezP.between=ezPlot(data=xdf, dv=.(abundance), between=.(treatment), sid = .(line_treat_id), do_lines=FALSE, x_lab='treatment', y_lab='abundance' , x=.(treatment))
#pdf(file = 'between_treatment_X22.pdf')
#print(ezP.between)
#dev.off()

#synergy
#ezP.between=ezPlot(data=xdf, dv=.(abundance), between=.(synergy), sid = .(line_treat_id), do_lines=FALSE, x_lab='synergy', y_lab='abundance' , x=.(synergy))
#pdf(file = 'between_synergy_X22.pdf')
#print(ezP.between)
#dev.off()

#Can also plot differences for "within" comparisons (i.e., time) if interested
ezP.within=ezPlot(data=xdf, dv=.(abundance), within=.(time_hrs), sid=.(line_treat_id), do_lines=TRUE, x_lab='time', y_lab='abundance' , x=.(time_hrs))
pdf(file = 'within_time_X22.pdf')
print(ezP.within)
dev.off()

#Plot Mixed
#time:treatment
ezP.mixed=ezPlot(data=xdf, dv=.(abundance), within=.(time_hrs), between=.(treatment), sid=.(line_treat_id), do_lines=FALSE, x_lab="time", y_lab="abundance" , x=.(time_hrs), split=.(treatment), split_lab = "treatment")
print(ezP.mixed)

#time:synergy
ezP.mixed2=ezPlot(data=xdf, dv=.(abundance), within=.(time_hrs), between=.(synergy), sid=.(line_treat_id), do_lines=FALSE, x_lab="time", y_lab="abundance" , x=.(time_hrs), split=.(synergy), split_lab = "Synergy")
print(ezP.mixed2)

#Need to figure out how to get legend/coloring to work better for example below.
#time:synergy:treatment
ezP.mixed2=ezPlot(data=xdf, dv=.(abundance), within=.(time_hrs), between=.(synergy,treatment), sid=.(line_treat_id), do_lines=FALSE, x_lab="time", y_lab="abundance" , x=.(time_hrs), split=.(synergy), split_lab = "Synergy")
print(ezP.mixed2)

