library("gplots")
library("heatmap.plus")
library(multtest)
library("ggplot2")

#Load drug annotations (target, type, DisplayName)
drug_annotations_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/CompoundsTargets.csv"
drug_annotations=read.csv(file=drug_annotations_file, row.names=1, as.is=c(1:4))

#Load clinical data
clin_data_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/ClinicalData 2 8 12-1.clean.txt"
rawclindata=read.table(clin_data_file,sep="\t",header=TRUE, row.names=1, as.is=c(1:20))

#Load subtype calls
subtype_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/UNC.TCGA.BRCA.ExpClusters.20110922_updatedSampleIDsV2.clean.txt"
subtypedata=read.table(subtype_file,sep="\t",header=TRUE, row.names=1, as.is=c(1:4))
#Exclude Normal tissue samples, and duplicate "-06A" samples
normals=which(subtypedata[,"Histology"]=="Normal Tissue")
subtypedata=subtypedata[-normals,]
duplicates=which(substr(rownames(subtypedata),13,16)=="-06A")
subtypedata=subtypedata[-duplicates,]
rownames(subtypedata)=substr(rownames(subtypedata),1,12) #Fix names to be consistent with clinical info file

#Load results file
path="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/Results/"
setwd(path)
#rawdata=read.table("TCGA_results_ExpMeth.txt",sep="\t",header=TRUE, as.is=c(1,2,3,4,6,12))
rawdata=read.table("TCGA_results_ExpMethCNV.txt",sep="\t",header=TRUE, as.is=c(1,2,3,4,6,12))

#Set output files
#bestdrugsfile="TCGA_results_ExpMeth_bestdrugs.txt"
bestdrugsfile="TCGA_results_ExpMethCNV_bestdrugs.txt"
#bestdrugsummaryfile="TCGA_results_ExpMeth_summary.txt"
bestdrugsummaryfile="TCGA_results_ExpMethCNV_summary.txt"
#heatmapfile="TCGA_results_ExpMeth_heatmap.pdf"
heatmapfile="TCGA_results_ExpMethCNV_heatmap.pdf"
#bestheatmapfile="TCGA_results_ExpMeth_heatmap_best.pdf"
bestheatmapfile="TCGA_results_ExpMethCNV_heatmap_best.pdf"
#statsummaryfile="TCGA_results_ExpMeth_statsummary.txt"
statsummaryfile="TCGA_results_ExpMethCNV_statsummary.txt"
patientdrugPercSensfile="TCGA_results_ExpMethCNV_PercSens.txt"
combinedfile="TCGA_results_ExpMethCNV_DrugProbs_w_ClinicalVars.txt"
LapVsHER2file="TCGA_results_ExpMethCNV_LapVsHER2.pdf"
TamoxVsERfile="TCGA_results_ExpMethCNV_TamoxVsER.pdf"
BIBW2992VsHER2file="TCGA_results_ExpMethCNV_BibwVsHER2.pdf"
probmatrixfile="TCGA_results_ExpMethCNV_prob_matrix.txt"

#Exclude duplicate "-06A" (metastatic) samples
duplicates=which(substr(rawdata[,"Sample"],13,16)=="-06A")
rawdata=rawdata[-duplicates,]

#Trim sample names to be consistent with clinical data
rawdata[,"Sample"]=substr(rawdata[,"Sample"],1,12)

#Filter out drugs that do not meet a specific modelAUC cutoff (used 0.6 as minimum in toolbox, but 0.7 might be better for summary figure and to be consistent with paper)
data=rawdata[which(rawdata[,"ModelAUC"]>=0.7),]

#Filter down to only patient with certain input data
#data=data[data[,"InputData"]=="Expression+Meth",]
data=data[data[,"InputData"]=="Expression+Meth+CNV",]

#Summarize most sensitive drug for all patients
#First, get list of patients and drugs 
patients=as.vector(unique(data[,"Sample"]))
drugs=unique(data[,"Drug.Compound"])

#Create dataframe to hold just best drug for each patient 
bestdrugdata = data.frame(cbind(patients, InputData=NA, Drug=NA, Model=NA, ModelAUC=NA, DataCombination=NA, NbModelVars=NA, PercModelVars=NA, WtPercModelVars=NA, Score=NA, Probability=NA, Sensitivity=NA), stringsAsFactors=FALSE, row.names=patients)

for (i in 1:length(patients)){
  patient=patients[i]
  patient_data=data[which(data[,"Sample"]==patient),]
  best_drug=patient_data[order(patient_data[,"Probability"], decreasing=TRUE)[1],]
  bestdrugdata[i,]=best_drug
}

#Write results to file
write.table(bestdrugdata, file=bestdrugsfile, sep="\t", row.names=FALSE)

#Summarize results
drug_counts=table(bestdrugdata[,"Drug"])
drug_names=names(drug_counts)
drug_percs=(drug_counts/sum(drug_counts))*100
drug_prob_mean=aggregate(as.numeric(bestdrugdata[,"Probability"]), by=list(bestdrugdata[,"Drug"]), mean)
drug_prob_sd=aggregate(as.numeric(bestdrugdata[,"Probability"]), by=list(bestdrugdata[,"Drug"]), sd)
bestdrugsummary=cbind(drug_names,as.vector(drug_counts),as.vector(drug_percs),drug_prob_mean[,"x"],drug_prob_sd[,"x"])
colnames(bestdrugsummary)=c("Drug compound","N samples","% samples","Mean predicted probability of response","Std predicted probability of response")
rank=order(as.numeric(bestdrugsummary[,"N samples"]), decreasing=TRUE)
bestdrugsummary=bestdrugsummary[rank,]

#Write summary to file
write.table(bestdrugsummary, file=bestdrugsummaryfile, sep="\t", row.names=FALSE)

#Create pie-chart to summarize
#pie(x=as.numeric(bestdrugsummary[,"N samples"]), labels=bestdrugsummary[,"Drug compound"])

#Create heatmap to summarize drugs vs patients
#First build matrix of probability values for drugs vs patients

#Create matrix to store drug vs patient data
prob_matrix=matrix(NA, nrow=length(patients), ncol=length(drugs), dimnames=list(patients,drugs))

for (i in 1:length(drugs)){
  drug=drugs[i]
  drug_data=data[which(data[,"Drug.Compound"]==drug),]
  prob_matrix[drug_data[,"Sample"],drug]=as.numeric(drug_data[,"Probability"])
}

#Determine columns (drugs) with no values (all missing) and remove
countNA=function(x){
length(which(is.na(x)))
}
missing_counts=apply(prob_matrix,2,countNA)
drugs_to_del=which(missing_counts==length(patients))
if(length(drugs_to_del)>0){
  prob_matrix=prob_matrix[,-drugs_to_del]
  drugs=colnames(prob_matrix)
}

#Write probability matrix to file
write.table(cbind(rownames(prob_matrix),prob_matrix), file=probmatrixfile, sep="\t", row.names=FALSE)

#Create Heatmaps
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}
### Heatmap 1 - All drugs that pass AUC cutoff and have some probability predictions
#Now with final set of patients and drugs, get associated annotations
subtypes=subtypedata[patients,"PAM50"]
drug_group=drug_annotations[drugs,4]
drug_names=drug_annotations[drugs,2]

#Set colors for color sidebars
drug_colors=drug_group
drug_colors[drug_colors==1]="red"
drug_colors[drug_colors==2]="blue"
drug_colors[drug_colors==3]="grey"

subtype_colors=subtypes
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="LumA"]="blue"
subtype_colors[subtype_colors=="LumB"]="cyan"
subtype_colors[subtype_colors=="Her2"]="pink"
subtype_colors[subtype_colors=="Claudin"]="yellow"
subtype_colors[subtype_colors=="Normal"]="green"

main_title="TCGA Breast Drug Sensitivity"
pdf(file=heatmapfile, height=7.5, width=10)
par(cex.main=1)
heatmap.2(t(prob_matrix), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(2,11), Rowv=TRUE, Colv=TRUE, ColSideColors=subtype_colors, RowSideColors=drug_colors, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=drug_names, cexRow=0.7, cexCol=0.9, col=rev(heat.colors(75)))
legend("topright",legend=c("Subtype","Basal","LumA","LumB","Her2","Claudin","Normal","","Drug type","Targeted","Chemo"), fill=c("white","red","blue","cyan","pink","yellow","green","white","white","red","blue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
#legend("topleft",legend=c("targeted","chemo"), fill=c("red","blue"), border=FALSE, bty="n", y.intersp = 0.5, cex=0.7)
dev.off()

### Heatmap 2 - All drugs that pass AUC cutoff and some minimal probability of response
#Exclude drugs with out at least some minimal probability of response in at least some patients
sens_thresh=0.65
Sens=which(apply(prob_matrix,2,max)>sens_thresh)
getPercSens=function(x){
 length(which(x>sens_thresh))/length(x)
}
getNumSens=function(x){
 length(which(x>sens_thresh))
}
NumSens=apply(prob_matrix,2,getNumSens)
PercSens=apply(prob_matrix,2,getPercSens)
prob_matrix2=prob_matrix[,Sens]
SensSummary065=cbind(NumSens,PercSens)
sens_thresh=0.50
NumSens=apply(prob_matrix,2,getNumSens)
PercSens=apply(prob_matrix,2,getPercSens)
SensSummary050=cbind(NumSens,PercSens)
write.table(cbind(SensSummary050,SensSummary065), file=patientdrugPercSensfile, sep="\t", row.names=TRUE)

#Exclude Bortezomib-batch2 to be consistent with Anneleen's tables
prob_matrix2=prob_matrix2[,-which(colnames(prob_matrix2)=="Bortezomib-batch2")]

patients=rownames(prob_matrix2)
drugs=colnames(prob_matrix2)

#Now with final set of patients and drugs, get associated annotations
subtypes=subtypedata[patients,"PAM50"]
ERstatus=rawclindata[patients,"ERStatus"]
PRstatus=rawclindata[patients,"PRStatus"]
HER2status=rawclindata[patients,"HER2FinalStatus"]
Tstatus=rawclindata[patients,"Tumor_T1_Coded"]
Nstatus=rawclindata[patients,"Node_Coded"]
Mstatus=rawclindata[patients,"Metastasis_Coded"]
drug_group=drug_annotations[drugs,4]
drug_names=drug_annotations[drugs,2]

#Set colors for color sidebars
drug_colors=drug_group
drug_colors[drug_colors==1]="red"
drug_colors[drug_colors==2]="blue"
drug_colors[drug_colors==3]="grey"

subtype_colors=subtypes
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="LumA"]="blue"
subtype_colors[subtype_colors=="LumB"]="cyan"
subtype_colors[subtype_colors=="Her2"]="pink"
subtype_colors[subtype_colors=="Claudin"]="yellow"
subtype_colors[subtype_colors=="Normal"]="green"

ER_colors=ERstatus
ER_colors[ER_colors=="Positive"]="black"
ER_colors[ER_colors=="Negative"]="white"
ER_colors[ER_colors!="black" & ER_colors!="white"]="grey"
ER_colors[is.na(ER_colors)]="grey"

PR_colors=PRstatus
PR_colors[PR_colors=="Positive"]="black"
PR_colors[PR_colors=="Negative"]="white"
PR_colors[PR_colors!="black" & PR_colors!="white"]="grey"
PR_colors[is.na(PR_colors)]="grey"

HER2_colors=HER2status
HER2_colors[HER2_colors=="Positive"]="black"
HER2_colors[HER2_colors=="Negative"]="white"
HER2_colors[HER2_colors!="black" & HER2_colors!="white"]="grey"
HER2_colors[is.na(HER2_colors)]="grey"

Tcolors=Tstatus
Tcolors[Tcolors=="T1"]="black"
Tcolors[Tcolors=="T_Other"]="white"
Tcolors[Tcolors!="black" & Tcolors!="white"]="grey"
Tcolors[is.na(Tcolors)]="grey"

Ncolors=Nstatus
Ncolors[Ncolors=="Positive"]="black"
Ncolors[Ncolors=="Negative"]="white"
Ncolors[Ncolors!="black" & Ncolors!="white"]="grey"
Ncolors[is.na(Ncolors)]="grey"

Mcolors=Mstatus
Mcolors[Mcolors=="Positive"]="black"
Mcolors[Mcolors=="Negative"]="white"
Mcolors[Mcolors!="black" & Mcolors!="white"]="grey"
Mcolors[is.na(Mcolors)]="grey"

#Use modified heatmap.2 command to allow multiple color side bars
source("C:/Users/Obi/Documents/My Dropbox/drug_predictors/Rscripts/heatmap.3.R")
rlab=cbind(drug_colors,drug_colors)
clab=cbind(subtype_colors,Mcolors,Ncolors,Tcolors,HER2_colors,PR_colors,ER_colors)
colnames(rlab)=c("","")
colnames(clab)=c("Subtype","M","N","T","HER2","PR","ER")

pdf(file=bestheatmapfile, height=7.5, width=10)
main_title="TCGA Breast Tumor Drug Response Predictions"
par(cex.main=1)
heatmap.3(t(prob_matrix2), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(2,14), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=drug_names, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=7, KeyValueName="Prob. Response",)
legend("topright",legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA"), fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
legend("left",legend=c("Targeted","Chemo"), fill=c("red","blue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#look for associations between predicted drug sensitivity and clinical variables (e.g., association between HER2 status and Lapatinib sensitivity)
patients=rownames(prob_matrix2)
drugs=colnames(prob_matrix2)
tests=c("ER","PR","HER2","T","N","M","Subtype","Stage")

#Create data objects to store results
stat_results=matrix(NA, nrow=length(drugs), ncol=length(tests), dimnames=list(drugs,tests))
stat_results_long=matrix(NA, nrow=length(drugs)*length(tests), ncol=3, dimnames=list( c(1:(length(drugs)*length(tests))), c("Drug","ClinVar","P-value") ) )
k=1
for (i in 1:length(drugs)){
 for (j in 1:length(tests)){
 stat_results_long[k,"Drug"]=drugs[i]
 stat_results_long[k,"ClinVar"]=tests[j]
 k=k+1
 }
}

for (i in 1:length(drugs)){
  drug=drugs[i]
  drug_data=prob_matrix2[patients,drug]
  ERstatus=rawclindata[patients,"ERStatus"]
  PRstatus=rawclindata[patients,"PRStatus"]
  HER2status=rawclindata[patients,"HER2FinalStatus"]
  Tstatus=rawclindata[patients,"Tumor_T1_Coded"]
  Nstatus=rawclindata[patients,"Node_Coded"]
  Mstatus=rawclindata[patients,"Metastasis_Coded"]
  subtypes=subtypedata[patients,"PAM50"]
  subtypedrugdata=as.data.frame(cbind(subtypes,as.numeric(drug_data)))
  stage=rawclindata[patients,"AJCC_Stage"]
  stage[which(stage=="Stage X")]=NA; stage[which(stage=="Not Available")]=NA #Clean up weird entries
  stagedata=as.data.frame(cbind(stage,as.numeric(drug_data)))

  ERstatus.wilcox=wilcox.test(x=drug_data[which(ERstatus=="Positive")], y=drug_data[which(ERstatus=="Negative")], alternative="two.sided")
  PRstatus.wilcox=wilcox.test(x=drug_data[which(PRstatus=="Positive")], y=drug_data[which(PRstatus=="Negative")], alternative="two.sided")
  HER2status.wilcox=wilcox.test(x=drug_data[which(HER2status=="Positive")], y=drug_data[which(HER2status=="Negative")], alternative="two.sided")
  Tstatus.wilcox=wilcox.test(x=drug_data[which(Tstatus=="T1")], y=drug_data[which(Tstatus=="T_Other")], alternative="two.sided")
  Nstatus.wilcox=wilcox.test(x=drug_data[which(Nstatus=="Positive")], y=drug_data[which(Nstatus=="Negative")], alternative="two.sided")
  Mstatus.wilcox=wilcox.test(x=drug_data[which(Mstatus=="Positive")], y=drug_data[which(Mstatus=="Negative")], alternative="two.sided")
  Subtype.aov=aov(drug_data~subtypes,data=subtypedrugdata)
  Stage.aov=aov(drug_data~stage,data=stagedata)

  stat_results[drug,"ER"]=ERstatus.wilcox$p.value
  stat_results[drug,"PR"]=PRstatus.wilcox$p.value
  stat_results[drug,"HER2"]=HER2status.wilcox$p.value
  stat_results[drug,"T"]=Tstatus.wilcox$p.value
  stat_results[drug,"N"]=Nstatus.wilcox$p.value
  stat_results[drug,"M"]=Mstatus.wilcox$p.value
  stat_results[drug,"Subtype"]=summary(Subtype.aov)[[1]][["Pr(>F)"]][1]
  stat_results[drug,"Stage"]=summary(Stage.aov)[[1]][["Pr(>F)"]][1]

  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="ER"),"P-value"]=ERstatus.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="PR"),"P-value"]=PRstatus.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="HER2"),"P-value"]=HER2status.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="T"),"P-value"]=Tstatus.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="N"),"P-value"]=Nstatus.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="M"),"P-value"]=Mstatus.wilcox$p.value
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="Subtype"),"P-value"]=summary(Subtype.aov)[[1]][["Pr(>F)"]][1]
  stat_results_long[which(stat_results_long[,"Drug"]==drug & stat_results_long[,"ClinVar"]=="Stage"),"P-value"]=summary(Stage.aov)[[1]][["Pr(>F)"]][1]
}

#Perform multiple testing corrections
#Correct p-values
pvalues=as.numeric(stat_results_long[,"P-value"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
stat_results_long=cbind(stat_results_long, pvalues_adj_orig_order[,2])
colnames(stat_results_long)[4]="BH"

write.table(stat_results_long, file=statsummaryfile, sep="\t", row.names=FALSE)

#Create matrix combining drug response probabilities with clinical data, write to file
combined_data=cbind(patients,prob_matrix2[patients,],rawclindata[patients,])
write.table(combined_data, file=combinedfile, sep="\t", row.names=FALSE)

#create plot illustrating probabilities of drug response for selected drugs versus clinical groups
#Tamoxifen vs ER status - limit to those with ER status reported
data=combined_data[which(combined_data[,"ERStatus"]=="Positive" | combined_data[,"ERStatus"]=="Negative"),c("Tamoxifen","ERStatus")]
pdf(file=TamoxVsERfile, height=7.5, width=10)
ggplot(data, aes(Tamoxifen, fill = ERStatus)) + geom_density(alpha = 0.2) + xlab("Tamoxifen Response Prob") + ylab("Density") + opts(title="Predicted Tamoxifen response by ER status")
dev.off()

#Lapatinib vs HER2 status
data=combined_data[which(combined_data[,"HER2FinalStatus"]=="Positive" | combined_data[,"HER2FinalStatus"]=="Negative"),c("GSK_Tykerb","HER2FinalStatus")]
pdf(file=LapVsHER2file, height=7.5, width=10)
ggplot(data, aes(GSK_Tykerb, fill = HER2FinalStatus)) + geom_density(alpha = 0.2) + xlab("Lapatinib Response Prob") + ylab("Density") + opts(title="Predicted Lapatinib response by HER2 status")
dev.off()

#BIBW2992 vs HER2 status
data=combined_data[which(combined_data[,"HER2FinalStatus"]=="Positive" | combined_data[,"HER2FinalStatus"]=="Negative"),c("BIBW2992","HER2FinalStatus")]
pdf(file=BIBW2992VsHER2file, height=7.5, width=10)
ggplot(data, aes(BIBW2992, fill = HER2FinalStatus)) + geom_density(alpha = 0.2) + xlab("BIBW2992 Response Prob") + ylab("Density") + opts(title="Predicted BIBW2992 response by HER2 status")
dev.off()

