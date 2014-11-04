#Load the appropriate libraries
library(affy)
library(gcrma)
library(multtest)
library(Biobase)
library(genefilter)
library("biomaRt")
library("gplots")
library("RColorBrewer")

#Set working directory for data files
setwd("/home/obig/Projects/Lymphoma_drugs/EZH2/CEL_FILES/DLBCL/SPEC_EZH2_GCB")
#setwd("/home/obig/Projects/Lymphoma_drugs/EZH2/CEL_FILES/DLBCL/SPECS")

Data=ReadAffy()

#Run GCRMA on Data to background correct, normalize and summarize expression data
gcrmaData=gcrma(Data)

#Set working data for results files
setwd("/home/obig/Projects/Lymphoma_drugs/EZH2/results/DLBCL")

#Write gcrma normalized data to file
#First, reduce number of decimal places
gcrma_values_formatted=format(exprs(gcrmaData), digits=5)
write.table(gcrma_values_formatted, file = "DLBCL_SPEC_EZH2_GCB_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "DLBCL_SPEC_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")

#Set working data for results files
setwd("/home/obig/Projects/Lymphoma_drugs/EZH2/results/DLBCL/diff_regulated")

#Get sample names and expression data from gcrma expression object
pheno=pData(gcrmaData)
X=exprs(gcrmaData)

#Create a vector of probe names for output later
var_names=rownames(X)

#Preliminary gene filtering might be a good idea. 
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation):
#At least 20% of samples should have raw intensity greater than 100
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_gcrmaData=gcrmaData[filt,]
filt_data=exprs(filt_gcrmaData)

#Create a vector of probe names after filtering for output later
var_names_filt=rownames(exprs(filt_gcrmaData))

#Then perform differential expression statistics and multiple testing correction
#Create a vector for the class labels
#EZH2 mt vs wt
class_labels=vector(length=35)
class_labels[c(1,2,4,5,13,15,19,28,29)]="mutant"
class_labels[c(3,6,7,8,9,10,11,12,14,16,17,18,20,21,22,23,24,25,26,27,30,31,32,33,34,35)]="wildtype"

#The Mann-Whitney test (or Wilcoxen Rank Test)="t.twosamp.equalvar" with "robust=TRUE"
#Also, note that B=1000 (the number of permutations) or higher is recommended.  But, for testing 100 is convenient.
MTP_results=MTP(X=filt_data, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP", robust=TRUE)

#Calculate basic summary stats (means, fold change, etc)
X=2^filt_data[,c(1,2,4,5,13,15,19,28,29)]
Y=2^filt_data[,c(3,6,7,8,9,10,11,12,14,16,17,18,20,21,22,23,24,25,26,27,30,31,32,33,34,35)]
class1_means=apply (X, 1, mean)
class2_means=apply (Y, 1, mean)
foldchanges=class2_means/class1_means
foldchanges[which(foldchanges<1)]=-1/foldchanges[which(foldchanges<1)]


#Output results.
MTP_summary=cbind(var_names_filt,class1_means,class2_means,foldchanges,MTP_results@statistic,MTP_results@estimate,MTP_results@rawp, MTP_results@adjp, MTP_results@reject)
colnames(MTP_summary)=c("probe", "mean_mt", "mean_wt", "FC", "statistic", "estimate", "rawp", "adjp", "reject")

#For only probes significant after multiple testing correction
#sigprobe_summary=MTP_summary[MTP_summary[,"reject"]=="TRUE",]

#For all 'significant' probes
sigprobe_summary=MTP_summary[MTP_summary[,"rawp"]<0.05,]

write.table(sigprobe_summary, file="DLBCL_SPEC_GCB_EZH2_mt_vs_wt.sigprobes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)


####################Annotate Probes to Genes or Proteins################################
#probe_ids=sigprobe_summary[,1]
probe_ids=MTP_summary[,1] #To map all filtered probes

#Create a mart object from ensembl Biomart
mart <- useMart("ensembl", "hsapiens_gene_ensembl")
attrbuts=listAttributes(mart)

########################################Uniprot mapping#########################################
#Get uniprot_swissprot_accessions for all probes
#annotations=getBM(attributes=c("affy_hg_u133_plus_2","uniprot_swissprot_accession"), filter="affy_hg_u133_plus_2", values=probe_ids, mart=mart, na.value="NA")   
#annotations=getBM(attributes=c("affy_hg_u133_plus_2","unified_uniprot_accession"), filter="affy_hg_u133_plus_2", values=probe_ids, mart=mart, na.value="NA")   
annotations=getBM(attributes=c("affy_hg_u133_plus_2","hgnc_symbol"), filter="affy_hg_u133_plus_2", values=probe_ids, mart=mart, na.value="NA")   

#Remove failed mappings (i.e. no uniprot found)
#annotations_notnull=annotations[which(annotations$uniprot_swissprot_accession!=""),]
annotations_notnull=annotations[which(annotations$hgnc_symbol!=""),]

#Remove redundant entries (where same probe is mapped to same uniprot - these are pointless)
annotations_notnull_unique=unique(annotations_notnull)

#Get number of mappings for each probe (we want only probes that map to one protein)
probe_counts_table=table(annotations_notnull_unique[,1]) 
probe_counts=probe_counts_table[annotations_notnull_unique[,1]]

#Add number of mappings per probe to annotations data
annotations_notnull_unique_mapcount=cbind(annotations_notnull_unique,probe_counts)

#If we order by probe_id we should see something sensible (i.e. that probes with two mappings to two uniprots have a mapcount of two)
#Check against live BioMart to make sure things look right
annotations_notnull_unique_mapcount[order(annotations_notnull_unique_mapcount[,1]),]

#Remove entries where one probe maps to multiple uniprot ids (still allows many probes to map to one uniprot)
annotations_notnull_unique_unambig=annotations_notnull_unique_mapcount[annotations_notnull_unique_mapcount[,3]==1,][,1:2]

#Write these mappings to a separate file.
#write.table(annotations_notnull_unique_unambig, file = "DLBCL_SPEC_GCB_EZH2_mt_vs_wt.sigprobes_rawp_2_Biomart_HGNC_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
write.table(annotations_notnull_unique_unambig, file = "DLBCL_SPEC_GCB_EZH2_mt_vs_wt.filtered_probes_2_Biomart_HGNC_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")


#Create a new summary file for just probe ids that also have hgnc gene names
sigprobes=rownames(MTP_summary[MTP_summary[,"rawp"]<0.05,])
sigprobes_w_hgnc=sigprobes[sigprobes%in%annotations_notnull_unique_unambig[,1]]
siggene_summary=sigprobe_summary[sigprobes_w_hgnc,]

#Look up HGNC ID for each probe
hgnc_ids=vector(length=length(sigprobes_w_hgnc))
for (i in 1:length(sigprobes_w_hgnc)){
  hgnc_ids[i]=annotations_notnull_unique_unambig[annotations_notnull_unique_unambig[,1]==sigprobes_w_hgnc[i],2]
}

siggene_summary=cbind(hgnc_ids,siggene_summary)
write.table(siggene_summary, file="DLBCL_SPEC_GCB_EZH2_mt_vs_wt.siggenes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)

###Heatmaps###
cond_colors=class_labels
cond_colors[cond_colors=="mutant"]="#bbbddc"
cond_colors[cond_colors=="wildtype"]="#756bb1"
sample_names=c("VAN_00042","VAN_00043","VAN_00045","VAN_00047","VAN_00048","VAN_00049","VAN_00050","VAN_00057","VAN_00064","VAN_00065","VAN_00133","VAN_00143","VAN_00144","VAN_00148","VAN_00150","VAN_00156","VAN_00567","VAN_00591","VAN_00594","VAN_00596","VAN_00598","VAN_00606","VAN_00608","VAN_00622","VAN_00633","VAN_00667","VAN_00675","VAN_00694","VAN_00697","VAN_00700","VAN_00702","VAN_00706","VAN_00793","VAN_01007","VAN_01008")

##only significant genes (before correction)
sigprobes=rownames(MTP_summary[MTP_summary[,"rawp"]<0.05,])
filt_data_sigprobes=filt_data[sigprobes,]
pdf("DLBCL_SPEC_GCB_EZH2_mt_vs_wt_sigprobes_heatmap.pdf")
mypalette<-brewer.pal(9,"Blues")
heatmap.2(filt_data_sigprobes, na.rm = TRUE, scale="none", symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=sample_names, col=mypalette, cexCol=0.75, ColSideColors=cond_colors, keysize=1.75)
dev.off()


##only significant genes (before correction) and with a hgnc gene name
sigprobes=rownames(MTP_summary[MTP_summary[,"rawp"]<0.05,])
sigprobes_w_hgnc=sigprobes[sigprobes%in%annotations_notnull_unique_unambig[,1]]
#Look up HGNC ID for each probe
hgnc_ids=vector(length=length(sigprobes_w_hgnc))
for (i in 1:length(sigprobes_w_hgnc)){
  hgnc_ids[i]=annotations_notnull_unique_unambig[annotations_notnull_unique_unambig[,1]==sigprobes_w_hgnc[i],2]
}
filt_data_sigprobes=filt_data[sigprobes_w_hgnc,]
pdf("DLBCL_SPEC_GCB_EZH2_mt_vs_wt_siggenes_heatmap.pdf")
mypalette<-brewer.pal(9,"Blues")
heatmap.2(filt_data_sigprobes, na.rm = TRUE, scale="none", symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=sample_names, col=mypalette, cexCol=0.75, ColSideColors=cond_colors)
dev.off()


##only significant genes (after correction)
sigprobes=rownames(MTP_summary[MTP_summary[,"reject"]=="TRUE",])
filt_data_sigprobes=filt_data[sigprobes,]
pdf("DLBCL_SPEC_GCB_EZH2_mt_vs_wt_sigprobes_corr_heatmap.pdf")
mypalette<-brewer.pal(9,"Blues")
heatmap.2(filt_data_sigprobes, na.rm = TRUE, scale="none", symkey=FALSE, density.info="none", trace="none", labRow=sigprobes, labCol=sample_names, col=mypalette, cexCol=0.75, ColSideColors=cond_colors)
dev.off()


##only significant genes (after correction) and with a hgnc gene name
sigprobes=rownames(MTP_summary[MTP_summary[,"reject"]=="TRUE",])
sigprobes_w_hgnc=sigprobes[sigprobes%in%annotations_notnull_unique_unambig[,1]]
#Look up HGNC ID for each probe
hgnc_ids=vector(length=length(sigprobes_w_hgnc))
for (i in 1:length(sigprobes_w_hgnc)){
  hgnc_ids[i]=annotations_notnull_unique_unambig[annotations_notnull_unique_unambig[,1]==sigprobes_w_hgnc[i],2]
}
filt_data_sigprobes=filt_data[sigprobes_w_hgnc,]
pdf("DLBCL_SPEC_GCB_EZH2_mt_vs_wt_siggenes_corr_heatmap.pdf")
mypalette<-brewer.pal(9,"Blues")
heatmap.2(filt_data_sigprobes, na.rm = TRUE, scale="none", symkey=FALSE, density.info="none", trace="none", labRow=hgnc_ids, labCol=sample_names, col=mypalette, cexCol=0.75, ColSideColors=cond_colors)
dev.off()
