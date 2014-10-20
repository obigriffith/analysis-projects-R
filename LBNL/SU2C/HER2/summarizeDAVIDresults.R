library("gplots")
library("heatmap.plus")

#Get all DAVID results
datadir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/DAVID_analysis/v2/PLvsLR/"
resultdir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/plotsTables/DAVID_analysis/overlap/"
resultfile="PLvsLR_DAVID_overlap_results.txt"
heatmapfile="PLvsLR_DAVID_overlap_results.pdf"



datafiles=list.files(datadir)
setwd(datadir)

#Create dataframe of correct format for first file, then add on all subsequent files
combined_data=read.table(datafiles[1], header=TRUE, sep="\t", as.is=c(1,2,6))
#Add column for rank and data comparison
combined_data=cbind(comparison=datafiles[1],rank=c(1:length(rownames(combined_data))),combined_data)
print(datafiles[1])

for (i in 2:length(datafiles)){
 data=read.table(datafiles[i], header=TRUE, sep="\t", as.is=c(1,2,6))
 data=cbind(comparison=datafiles[i],rank=c(1:length(rownames(data))),data)
 combined_data=rbind(combined_data,data)
 print(datafiles[i])
}

setwd(resultdir)

#Clean up comparisons names (expecting file names like: PLvsLR_HCC1954_DAVID.txt)
combined_data[,"comparison"]=gsub("PLvsLR_","",combined_data[,"comparison"])
combined_data[,"comparison"]=gsub("_DAVID.txt","",combined_data[,"comparison"])

#Eliminate duplicated terms (e.g., "disulfide bond") within model. If repeated Term from different categories keeps most significant one.
combined_data=combined_data[-which(duplicated(paste(combined_data[,"comparison"], combined_data[,"Term"], sep="____"))),]

#Summarize results term by term
term_counts=table(combined_data[,"Term"]) #number of datasets with each term

#Calculate ranksum, mean p-value, etc
term_average_pvalue=aggregate.data.frame(combined_data[,"PValue"], by=list(combined_data[,"Term"]), mean)
term_average_FE=aggregate.data.frame(combined_data[,"Fold.Enrichment"], by=list(combined_data[,"Term"]), mean)
term_average_Benjamini=aggregate.data.frame(combined_data[,"Benjamini"], by=list(combined_data[,"Term"]), mean)
term_average_FDR=aggregate.data.frame(combined_data[,"FDR"], by=list(combined_data[,"Term"]), mean)
term_ranksum=aggregate.data.frame(combined_data[,"rank"], by=list(combined_data[,"Term"]), sum)

#Term to Category mapping
TermToCategory=matrix(unlist(strsplit(unique(paste(combined_data[,"Category"], combined_data[,"Term"], sep="____")), split="____")), ncol=2, byrow=TRUE)
rownames(TermToCategory)=TermToCategory[,2]

#Create shorter versions of Term names
TermToCategory=cbind(TermToCategory,TermToCategory[,2])
GOTermsind=which(TermToCategory[,1]=="GOTERM_BP_FAT" | TermToCategory[,1]=="GOTERM_CC_FAT" | TermToCategory[,1]=="GOTERM_MF_FAT")
TermToCategory[GOTermsind,3]=substring(TermToCategory[GOTermsind,2],first=12)
IPRTermsind=which(TermToCategory[,1]=="INTERPRO")
TermToCategory[IPRTermsind,3]=substring(TermToCategory[IPRTermsind,2],first=11)
SMARTTermsind=which(TermToCategory[,1]=="SMART")
TermToCategory[SMARTTermsind,3]=substring(TermToCategory[SMARTTermsind,2],first=9)
PIRTermsind=which(TermToCategory[,1]=="PIR_SUPERFAMILY")
TermToCategory[PIRTermsind,3]=substring(TermToCategory[PIRTermsind,2],first=13)
KEGGTermsind=which(TermToCategory[,1]=="KEGG_PATHWAY")
TermToCategory[KEGGTermsind,3]=substring(TermToCategory[KEGGTermsind,2],first=10)
BBIDTermsind=which(TermToCategory[,1]=="BBID")
TermToCategory[BBIDTermsind,3]=substring(TermToCategory[BBIDTermsind,2],first=4)

#For each term, create combined list of comparisons
Term_complist=aggregate.data.frame(combined_data[,"comparison"], by=list(combined_data[,"Term"]), FUN=paste, collapse=",")

#For each term create combined gene list, collapsed down to unique genes
Term_genelist=aggregate.data.frame(combined_data[,"Genes"], by=list(combined_data[,"Term"]), FUN=paste, collapse=", ")
for(i in 1:length(Term_genelist[,2])){
 Term_genelist[i,2]=paste(sort(unique(strsplit(Term_genelist[i,2], split=", ")[[1]])), collapse=",")
}

#Write results to file
results=cbind(names(term_counts),TermToCategory[names(term_counts),1],term_counts,Term_complist[,"x"],term_average_pvalue[,"x"],term_average_FE[,"x"],term_average_Benjamini[,"x"],term_average_FDR[,"x"],term_ranksum[,"x"],Term_genelist[,"x"])
colnames(results)=c("Term","Category","ModelCount","Models","avg_p","avg_FE","avg_Benj","avg_FDR","ranksum","Genes")
#reorder based on avg FDR
results=results[order(as.numeric(results[,"avg_FDR"])),]
write.table(results, file=resultfile, quote=FALSE, sep="\t", row.names=FALSE)

#Limit to just results with an average FDR < 1
results_filt=results[which(as.numeric(results[,"avg_FDR"])<5),]
results_display=results_filt[,c("Term","Category","Models","avg_FDR")]

#Create matrix of average FDR values for Term vs Cell line
terms=unique(combined_data[,"Term"])
models=unique(combined_data[,"comparison"])
FDRdata=matrix(NA, nrow=length(terms), ncol=length(models), byrow=FALSE, dimnames=list(terms,models))
for (i in 1:length(rownames(combined_data))){
 FDRdata[combined_data[i,"Term"],combined_data[i,"comparison"]]=combined_data[i,"FDR"]
}

#Create figure to summarize the results
#First cut down to just filtered set of terms
terms_filt=rownames(results_filt)
FDRdata_filt=FDRdata[terms_filt,]
term_names_filt=TermToCategory[terms_filt,3]
term_categories_filt=TermToCategory[terms_filt,1]

#Set NA values for FDR to 100
FDRdata_filt[is.na(FDRdata_filt)]=100

#models: "AU565", "BT474AZ", "HCC1954", "HCC202", "MD361", "SKBR3", "UACC812"
ERstatus=c("0","1","0","0","1","0","1")
PRstatus=c("0","1","0","0","0","0","0")
EGFRstatus=c("3","1","3","1","1","3","1")
subtype=c("Luminal","Luminal","Basal","Luminal","Luminal","Luminal","Luminal")
sampledata=data.frame(ER=ERstatus,PR=PRstatus,EGFR=EGFRstatus,subtype=subtype, row.names=models)

#Re-order data and sample data for arbitrary sample order
newsampleorder=c("HCC1954","AU565","SKBR3","HCC202","UACC812","MD361","BT474AZ")
FDRdata_filt=FDRdata_filt[,newsampleorder]
sampledata=sampledata[newsampleorder,]

subtype_colors=as.vector(sampledata[,"subtype"])
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Luminal"]="cyan"

ER_colors=as.vector(sampledata[,"ER"])
ER_colors[ER_colors=="0"]="#FFFFFF"
ER_colors[ER_colors=="1"]="#3182BD"

PR_colors=as.vector(sampledata[,"PR"])
PR_colors[PR_colors=="0"]="#FFFFFF"
PR_colors[PR_colors=="1"]="#3182BD"

EGFR_colors=as.vector(sampledata[,"EGFR"])
EGFR_colors[EGFR_colors=="0"]="#FFFFFF"
EGFR_colors[EGFR_colors=="1"]="#DEEBF7"
EGFR_colors[EGFR_colors=="2"]="#9ECAE1"
EGFR_colors[EGFR_colors=="3"]="#3182BD"

categories=unique(term_categories_filt)
category_numbers=c(1:length(categories))
names(category_numbers)=categories
category_colors=rainbow(length(categories))
term_category_colors=category_colors[category_numbers[term_categories_filt]]

rlab=cbind(term_category_colors,term_category_colors)
clab=cbind(subtype_colors,EGFR_colors,ER_colors,PR_colors)
colnames(rlab)=c("Category","")
colnames(clab)=c("Subtype","EGFR","ER","PR")

#Use modified heatmap.2 command to allow multiple color side bars
source("C:/Users/Obi/Documents/My Dropbox/drug_predictors/Rscripts/heatmap.3.R")
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Transform data to -log10 for easier visualization. Divide by 100 to counteract FDR values having been turned into percents by DAVID
FDRdata_filt_trans = -log10(FDRdata_filt/100)

pdf(file=heatmapfile, height=14, width=8.5)
main_title="Signif. Pathways for HER2 models"
par(cex.main=1)
breaks=c(0,0.5,25,50,75,100)
colors=c("white","#807DBA","#6A51A3","#54278F","#3F007D")
heatmap.3(FDRdata_filt_trans, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="none", margins=c(6,24), Rowv=FALSE, Colv=FALSE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.7, cexCol=1, breaks=breaks, col=colors, labRow=term_names_filt, NumColSideColors=3, NumRowSideColors=2, KeyValueName="-log10(FDR)")
legend("topright",legend=c("Basal","Luminal","","Positive/High","Low","Negative"), fill=c("red","cyan","white","#3182BD","#DEEBF7","white"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
legend("left",legend=categories, fill=category_colors, border=FALSE, bty="n", y.intersp = 0.8, cex=0.6)
dev.off()

#Make binary version of data: 0 not sig, 1 sig
FDRdata_filt_trans_bin=FDRdata_filt_trans
FDRdata_filt_trans_bin[FDRdata_filt_trans_bin>0]=1
sort(apply(FDRdata_filt_trans_bin[,5:7],1,sum)-apply(FDRdata_filt_trans_bin[,1:4],1,sum))

