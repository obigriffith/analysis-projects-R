#Get all sigDE csv files
#datadir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/sigDE/"
#datadir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/Her2_sigDe/"
datadir="/Users/obigriffith/Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/Her2_sigDe/"
#overlapdir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/overlap_analysis/"
overlapdir="/Users/obigriffith/Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/overlap_analysis/"

datafiles=list.files(datadir)
setwd(datadir)

#Create dataframe of correct format for first file, then add on all subsequent files
combined_data=read.table(datafiles[1], header=TRUE)
for (i in 2:length(datafiles)){
 print(datafiles[i])
 data=read.table(datafiles[i], header=TRUE)
 combined_data=rbind(combined_data,data)
}

setwd(overlapdir)

##Find overlap between all 5 models for P vs LR (or LLR)##
#Get sample rows for correct comparisons
BT474AZ_PvsLLR = which(combined_data[,"sample_1"]=="BT474AZ_P" & combined_data[,"sample_2"]=="BT474AZ_LLR")
HCC1954_PvsLR = which(combined_data[,"sample_1"]=="HCC1954_P" & combined_data[,"sample_2"]=="HCC1954_LR")
SKBR3_PvsLR = which(combined_data[,"sample_1"]=="SKBR3_P" & combined_data[,"sample_2"]=="SKBR3_LR")
UACC812_PvsLR = which(combined_data[,"sample_1"]=="UACC812_P" & combined_data[,"sample_2"]=="UACC812_LR")
HCC202_PvsLR = which(combined_data[,"sample_1"]=="HCC202_P" & combined_data[,"sample_2"]=="HCC202_LR")

#Extract unique gene symbols
BT474AZ_PvsLLR_genes=unique(as.vector(combined_data[BT474AZ_PvsLLR,"hg19.kgXref.geneSymbol"]))
HCC1954_PvsLR_genes=unique(as.vector(combined_data[HCC1954_PvsLR,"hg19.kgXref.geneSymbol"]))
SKBR3_PvsLR_genes=unique(as.vector(combined_data[SKBR3_PvsLR,"hg19.kgXref.geneSymbol"]))
UACC812_PvsLR_genes=unique(as.vector(combined_data[UACC812_PvsLR,"hg19.kgXref.geneSymbol"]))
HCC202_PvsLR_genes=unique(as.vector(combined_data[HCC202_PvsLR,"hg19.kgXref.geneSymbol"]))
All_PvsLR_genes=c(BT474AZ_PvsLLR_genes,HCC1954_PvsLR_genes,SKBR3_PvsLR_genes,UACC812_PvsLR_genes,HCC202_PvsLR_genes)

#Combine into a single file with format suitable for Venn diagram
BT474AZ_PvsLLR_venn=cbind(BT474AZ_PvsLLR_genes,rep("BT474AZ_PvsLLR",length(BT474AZ_PvsLLR_genes)))
HCC1954_PvsLR_venn=cbind(HCC1954_PvsLR_genes,rep("HCC1954_PvsLR",length(HCC1954_PvsLR_genes)))
SKBR3_PvsLR_venn=cbind(SKBR3_PvsLR_genes,rep("SKBR3_PvsLR",length(SKBR3_PvsLR_genes)))
UACC812_PvsLR_venn=cbind(UACC812_PvsLR_genes,rep("UACC812_PvsLR",length(UACC812_PvsLR_genes)))
HCC202_PvsLR_venn=cbind(HCC202_PvsLR_genes,rep("HCC202_PvsLR",length(HCC202_PvsLR_genes)))

All_PvsLR_venn=rbind(BT474AZ_PvsLLR_venn,HCC1954_PvsLR_venn,SKBR3_PvsLR_venn,UACC812_PvsLR_venn, HCC202_PvsLR_venn)
write.table(All_PvsLR_venn, file="All_PvsLR_venn_list.txt", sep="\t", row.names=FALSE, col.names=FALSE)

All_PvsLR_overlap_counts=sort(table(All_PvsLR_venn[,1]))
write.table(All_PvsLR_overlap_counts, file="All_PvsLR_overlap_counts.txt", sep="\t", row.names=TRUE, col.names=FALSE)


##Find overlap between all 5 models for P vs TR)##
#Get sample rows for correct comparisons
BT474AZ_PvsTR = which(combined_data[,"sample_1"]=="BT474AZ_P" & combined_data[,"sample_2"]=="BT474AZ_TR")
HCC1954_PvsTR = which(combined_data[,"sample_1"]=="HCC1954_P" & combined_data[,"sample_2"]=="HCC1954_TR")
SKBR3_PvsTR = which(combined_data[,"sample_1"]=="SKBR3_P" & combined_data[,"sample_2"]=="SKBR3_TR")
UACC812_PvsTR = which(combined_data[,"sample_1"]=="UACC812_P" & combined_data[,"sample_2"]=="UACC812_TR")
HCC202_PvsTR = which(combined_data[,"sample_1"]=="HCC202_P" & combined_data[,"sample_2"]=="HCC202_TR")

#Extract unique gene symbols
BT474AZ_PvsTR_genes=unique(as.vector(combined_data[BT474AZ_PvsTR,"hg19.kgXref.geneSymbol"]))
HCC1954_PvsTR_genes=unique(as.vector(combined_data[HCC1954_PvsTR,"hg19.kgXref.geneSymbol"]))
SKBR3_PvsTR_genes=unique(as.vector(combined_data[SKBR3_PvsTR,"hg19.kgXref.geneSymbol"]))
UACC812_PvsTR_genes=unique(as.vector(combined_data[UACC812_PvsTR,"hg19.kgXref.geneSymbol"]))
HCC202_PvsTR_genes=unique(as.vector(combined_data[HCC202_PvsTR,"hg19.kgXref.geneSymbol"]))
All_PvsTR_genes=c(BT474AZ_PvsTR_genes,HCC1954_PvsTR_genes,SKBR3_PvsTR_genes,UACC812_PvsTR_genes,HCC202_PvsTR_genes)

#Combine into a single file with format suitable for Venn diagram
BT474AZ_PvsTR_venn=cbind(BT474AZ_PvsTR_genes,rep("BT474AZ_PvsTR",length(BT474AZ_PvsTR_genes)))
HCC1954_PvsTR_venn=cbind(HCC1954_PvsTR_genes,rep("HCC1954_PvsTR",length(HCC1954_PvsTR_genes)))
SKBR3_PvsTR_venn=cbind(SKBR3_PvsTR_genes,rep("SKBR3_PvsTR",length(SKBR3_PvsTR_genes)))
UACC812_PvsTR_venn=cbind(UACC812_PvsTR_genes,rep("UACC812_PvsTR",length(UACC812_PvsTR_genes)))
HCC202_PvsTR_venn=cbind(HCC202_PvsTR_genes,rep("HCC202_PvsTR",length(HCC202_PvsTR_genes)))

All_PvsTR_venn=rbind(BT474AZ_PvsTR_venn,HCC1954_PvsTR_venn,SKBR3_PvsTR_venn,UACC812_PvsTR_venn,HCC202_PvsTR_venn)
write.table(All_PvsTR_venn, file="All_PvsTR_venn_list.txt", sep="\t", row.names=FALSE, col.names=FALSE)

All_PvsTR_overlap_counts=sort(table(All_PvsTR_venn[,1]))
write.table(All_PvsTR_overlap_counts, file="All_PvsTR_overlap_counts.txt", sep="\t", row.names=TRUE, col.names=FALSE)


##Find overlap between all 5 models for P vs LTR)##
#Get sample rows for correct comparisons
BT474AZ_PvsLTR = which(combined_data[,"sample_1"]=="BT474AZ_P" & combined_data[,"sample_2"]=="BT474AZ_LTR")
HCC1954_PvsLTR = which(combined_data[,"sample_1"]=="HCC1954_P" & combined_data[,"sample_2"]=="HCC1954_LTR")
SKBR3_PvsLTR = which(combined_data[,"sample_1"]=="SKBR3_P" & combined_data[,"sample_2"]=="SKBR3_LTR")
UACC812_PvsLTR = which(combined_data[,"sample_1"]=="UACC812_P" & combined_data[,"sample_2"]=="UACC812_LTR")
HCC202_PvsLTR = which(combined_data[,"sample_1"]=="HCC202_P" & combined_data[,"sample_2"]=="HCC202_LTR")

#Extract unique gene symbols
BT474AZ_PvsLTR_genes=unique(as.vector(combined_data[BT474AZ_PvsLTR,"hg19.kgXref.geneSymbol"]))
HCC1954_PvsLTR_genes=unique(as.vector(combined_data[HCC1954_PvsLTR,"hg19.kgXref.geneSymbol"]))
SKBR3_PvsLTR_genes=unique(as.vector(combined_data[SKBR3_PvsLTR,"hg19.kgXref.geneSymbol"]))
UACC812_PvsLTR_genes=unique(as.vector(combined_data[UACC812_PvsLTR,"hg19.kgXref.geneSymbol"]))
HCC202_PvsLTR_genes=unique(as.vector(combined_data[HCC202_PvsLTR,"hg19.kgXref.geneSymbol"]))
All_PvsLTR_genes=c(BT474AZ_PvsLTR_genes,HCC1954_PvsLTR_genes,SKBR3_PvsLTR_genes,UACC812_PvsLTR_genes,HCC202_PvsLTR_genes)

#Combine into a single file with format suitable for Venn diagram
BT474AZ_PvsLTR_venn=cbind(BT474AZ_PvsLTR_genes,rep("BT474AZ_PvsLTR",length(BT474AZ_PvsLTR_genes)))
HCC1954_PvsLTR_venn=cbind(HCC1954_PvsLTR_genes,rep("HCC1954_PvsLTR",length(HCC1954_PvsLTR_genes)))
SKBR3_PvsLTR_venn=cbind(SKBR3_PvsLTR_genes,rep("SKBR3_PvsLTR",length(SKBR3_PvsLTR_genes)))
UACC812_PvsLTR_venn=cbind(UACC812_PvsLTR_genes,rep("UACC812_PvsLTR",length(UACC812_PvsLTR_genes)))
HCC202_PvsLTR_venn=cbind(HCC202_PvsLTR_genes,rep("HCC202_PvsLTR",length(HCC202_PvsLTR_genes)))

All_PvsLTR_venn=rbind(BT474AZ_PvsLTR_venn,HCC1954_PvsLTR_venn,SKBR3_PvsLTR_venn,UACC812_PvsLTR_venn,HCC202_PvsLTR_venn)
write.table(All_PvsLTR_venn, file="All_PvsLTR_venn_list.txt", sep="\t", row.names=FALSE, col.names=FALSE)

All_PvsLTR_overlap_counts=sort(table(All_PvsLTR_venn[,1]))
write.table(All_PvsLTR_overlap_counts, file="All_PvsLTR_overlap_counts.txt", sep="\t", row.names=TRUE, col.names=FALSE)


##Find overlap between all 5 models for TR vs LR)##
#Get sample rows for correct comparisons
BT474AZ_TRvsLLR = which(combined_data[,"sample_1"]=="BT474AZ_TR" & combined_data[,"sample_2"]=="BT474AZ_LLR")
HCC1954_TRvsLR = which(combined_data[,"sample_1"]=="HCC1954_TR" & combined_data[,"sample_2"]=="HCC1954_LR")
SKBR3_TRvsLR = which(combined_data[,"sample_1"]=="SKBR3_TR" & combined_data[,"sample_2"]=="SKBR3_LR")
UACC812_TRvsLR = which(combined_data[,"sample_1"]=="UACC812_TR" & combined_data[,"sample_2"]=="UACC812_LR")
HCC202_TRvsLR = which(combined_data[,"sample_1"]=="HCC202_TR" & combined_data[,"sample_2"]=="HCC202_LR")

#Extract unique gene symbols
BT474AZ_TRvsLLR_genes=unique(as.vector(combined_data[BT474AZ_TRvsLLR,"hg19.kgXref.geneSymbol"]))
HCC1954_TRvsLR_genes=unique(as.vector(combined_data[HCC1954_TRvsLR,"hg19.kgXref.geneSymbol"]))
SKBR3_TRvsLR_genes=unique(as.vector(combined_data[SKBR3_TRvsLR,"hg19.kgXref.geneSymbol"]))
UACC812_TRvsLR_genes=unique(as.vector(combined_data[UACC812_TRvsLR,"hg19.kgXref.geneSymbol"]))
HCC202_TRvsLR_genes=unique(as.vector(combined_data[HCC202_TRvsLR,"hg19.kgXref.geneSymbol"]))
All_TRvsLR_genes=c(BT474AZ_TRvsLLR_genes,HCC1954_TRvsLR_genes,SKBR3_TRvsLR_genes,UACC812_TRvsLR_genes,HCC202_TRvsLR_genes)

#Combine into a single file with format suitable for Venn diagram
BT474AZ_TRvsLLR_venn=cbind(BT474AZ_TRvsLLR_genes,rep("BT474AZ_TRvsLLR",length(BT474AZ_TRvsLLR_genes)))
HCC1954_TRvsLR_venn=cbind(HCC1954_TRvsLR_genes,rep("HCC1954_TRvsLR",length(HCC1954_TRvsLR_genes)))
SKBR3_TRvsLR_venn=cbind(SKBR3_TRvsLR_genes,rep("SKBR3_TRvsLR",length(SKBR3_TRvsLR_genes)))
UACC812_TRvsLR_venn=cbind(UACC812_TRvsLR_genes,rep("UACC812_TRvsLR",length(UACC812_TRvsLR_genes)))
HCC202_TRvsLR_venn=cbind(HCC202_TRvsLR_genes,rep("HCC202_TRvsLR",length(HCC202_TRvsLR_genes)))

All_TRvsLR_venn=rbind(BT474AZ_TRvsLLR_venn,HCC1954_TRvsLR_venn,SKBR3_TRvsLR_venn,UACC812_TRvsLR_venn,HCC202_TRvsLR_venn)
write.table(All_TRvsLR_venn, file="All_TRvsLR_venn_list.txt", sep="\t", row.names=FALSE, col.names=FALSE)

All_TRvsLR_overlap_counts=sort(table(All_TRvsLR_venn[,1]))
write.table(All_TRvsLR_overlap_counts, file="All_TRvsLR_overlap_counts.txt", sep="\t", row.names=TRUE, col.names=FALSE)






##Find overlap between all 5 models for TR vs LTR)##
#Get sample rows for correct comparisons
BT474AZ_TRvsLTR = which(combined_data[,"sample_1"]=="BT474AZ_TR" & combined_data[,"sample_2"]=="BT474AZ_LTR")
HCC1954_TRvsLTR = which(combined_data[,"sample_1"]=="HCC1954_TR" & combined_data[,"sample_2"]=="HCC1954_LTR")
SKBR3_TRvsLTR = which(combined_data[,"sample_1"]=="SKBR3_TR" & combined_data[,"sample_2"]=="SKBR3_LTR")
UACC812_TRvsLTR = which(combined_data[,"sample_1"]=="UACC812_TR" & combined_data[,"sample_2"]=="UACC812_LTR")
HCC202_TRvsLTR = which(combined_data[,"sample_1"]=="HCC202_TR" & combined_data[,"sample_2"]=="HCC202_LTR")

#Extract unique gene symbols
BT474AZ_TRvsLTR_genes=unique(as.vector(combined_data[BT474AZ_TRvsLTR,"hg19.kgXref.geneSymbol"]))
HCC1954_TRvsLTR_genes=unique(as.vector(combined_data[HCC1954_TRvsLTR,"hg19.kgXref.geneSymbol"]))
SKBR3_TRvsLTR_genes=unique(as.vector(combined_data[SKBR3_TRvsLTR,"hg19.kgXref.geneSymbol"]))
UACC812_TRvsLTR_genes=unique(as.vector(combined_data[UACC812_TRvsLTR,"hg19.kgXref.geneSymbol"]))
HCC202_TRvsLTR_genes=unique(as.vector(combined_data[HCC202_TRvsLTR,"hg19.kgXref.geneSymbol"]))
All_TRvsLTR_genes=c(BT474AZ_TRvsLTR_genes,HCC1954_TRvsLTR_genes,SKBR3_TRvsLTR_genes,UACC812_TRvsLTR_genes,HCC202_TRvsLTR_genes)

#Combine into a single file with format suitable for Venn diagram
BT474AZ_TRvsLTR_venn=cbind(BT474AZ_TRvsLTR_genes,rep("BT474AZ_TRvsLTR",length(BT474AZ_TRvsLTR_genes)))
HCC1954_TRvsLTR_venn=cbind(HCC1954_TRvsLTR_genes,rep("HCC1954_TRvsLTR",length(HCC1954_TRvsLTR_genes)))
SKBR3_TRvsLTR_venn=cbind(SKBR3_TRvsLTR_genes,rep("SKBR3_TRvsLTR",length(SKBR3_TRvsLTR_genes)))
UACC812_TRvsLTR_venn=cbind(UACC812_TRvsLTR_genes,rep("UACC812_TRvsLTR",length(UACC812_TRvsLTR_genes)))
HCC202_TRvsLTR_venn=cbind(HCC202_TRvsLTR_genes,rep("HCC202_TRvsLTR",length(HCC202_TRvsLTR_genes)))

All_TRvsLTR_venn=rbind(BT474AZ_TRvsLTR_venn,HCC1954_TRvsLTR_venn,SKBR3_TRvsLTR_venn,UACC812_TRvsLTR_venn,HCC202_TRvsLTR_venn)
write.table(All_TRvsLTR_venn, file="All_TRvsLTR_venn_list.txt", sep="\t", row.names=FALSE, col.names=FALSE)

All_TRvsLTR_overlap_counts=sort(table(All_TRvsLTR_venn[,1]))
write.table(All_TRvsLTR_overlap_counts, file="All_TRvsLTR_overlap_counts.txt", sep="\t", row.names=TRUE, col.names=FALSE)





