#Load the appropriate libraries
library(multtest)
library(genefilter)
library("biomaRt")
library("gplots")
library("RColorBrewer")

#Load expression data
setwd("/home/obig/Projects/Colon_Cancer_Resistance/exon_array_data")
datafile="EnsEMBL_Genes_v53_iterplier_backgroundCorrected_sketchQuantiles.summary.clean.txt"
rawdata=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", row.names=2)

#Put aside probe names
alexa_ids=as.data.frame(rawdata[,1])
data=rawdata[,2:length(colnames(rawdata))]
rownames(alexa_ids)=rownames(rawdata)

#Create a vector of experiment and probe/gene names
all_exp_names=colnames(data)
all_gene_names=rownames(data)

#Preliminary gene filtering might be a good idea.
#Take values and filter out any genes according to following criteria (recommended in multtest/MTP documentation):
#At least 20% of samples should have raw intensity greater than 100
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(data,ffun)
filt_data=data[filt,]

#Create index_vectors for each cell_line/condition
MIP101=c(1,3,5,7)
MIP5FUR=c(2,4,6,8)
RKO=c(9,11,12,13)
RKO5FUR=c(10,14,15,16)
HCT=c(17,18,19)
HCT5FUR=c(20,21,22)
MIP101_5FU=c(23,24,25)
MIP5FUR_5FU=c(26,27,28)

#Perform each comparison
###Choose only one of the following comparisons before any subsequent commands###
#MIP101_vs_MIP5FUR
setwd("/home/obig/Projects/Colon_Cancer_Resistance/analysis/MIP101_vs_MIP5FUR")
class1="MIP101"
class2="MIP5FUR"
filt_data_sub=filt_data[,sort(c(MIP101,MIP5FUR))]
X=filt_data[,MIP101]
Y=filt_data[,MIP5FUR]
gene_names_filt=rownames(filt_data_sub)
exp_names_filt=colnames(filt_data_sub)
class_labels=vector(length=8)
class_labels[MIP101]="MIP101"
class_labels[MIP5FUR]="MIP5FUR"

#RKO vs RKO5FUR
setwd("/home/obig/Projects/Colon_Cancer_Resistance/analysis/RKO_vs_RKO5FUR")
class1="RKO"
class2="RKO5FUR"
filt_data_sub=filt_data[,sort(c(RKO,RKO5FUR))]
X=filt_data[,RKO]
Y=filt_data[,RKO5FUR]
gene_names_filt=rownames(filt_data_sub)
exp_names_filt=colnames(filt_data_sub)
class_labels=vector(length=8)
class_labels[c(1,3,4,5)]="RKO"
class_labels[c(2,6,7,8)]="RKO5FUR"

#HCT vs HCT5FUR
setwd("/home/obig/Projects/Colon_Cancer_Resistance/analysis/HCT_vs_HCT5FUR")
class1="HCT"
class2="HCT5FUR"
filt_data_sub=filt_data[,sort(c(HCT,HCT5FUR))]
X=filt_data[,HCT]
Y=filt_data[,HCT5FUR]
gene_names_filt=rownames(filt_data_sub)
exp_names_filt=colnames(filt_data_sub)
class_labels=vector(length=6)
class_labels[c(1,2,3)]="HCT"
class_labels[c(4,5,6)]="HCT5FUR"

#MIP101_5FU vs MIP5FUR_5FU
setwd("/home/obig/Projects/Colon_Cancer_Resistance/analysis/MIP101_5FU_vs_MIP5FUR_5FU")
class1="MIP101_5FU"
class2="MIP5FUR_5FU"
filt_data_sub=filt_data[,sort(c(MIP101_5FU,MIP5FUR_5FU))]
X=filt_data[,MIP101_5FU]
Y=filt_data[,MIP5FUR_5FU]
gene_names_filt=rownames(filt_data_sub)
exp_names_filt=colnames(filt_data_sub)
class_labels=vector(length=6)
class_labels[c(1,2,3)]="MIP101_5FU"
class_labels[c(4,5,6)]="MIP5FUR_5FU"


##################################################################################

#Once filtered down to final set, retrieve probe ids for remaining ENSG ids
alexa_ids_filt=vector(length=length(gene_names_filt))
for (i in 1:length(gene_names_filt)){
alexa_ids_filt[i]=alexa_ids[gene_names_filt[i],]
}

#Calculate basic summary stats (means, fold change, etc)
class1_means=apply (X, 1, mean)
class2_means=apply (Y, 1, mean)
foldchanges=class2_means/class1_means
foldchanges[which(foldchanges<1)]=-1/foldchanges[which(foldchanges<1)]

#Then perform differential expression statistics and multiple testing correction for each comparison
MTP_results=MTP(X=filt_data_sub, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP")

#Output results.
MTP_summary=cbind(gene_names_filt,alexa_ids_filt,class1_means,class2_means,foldchanges,MTP_results@rawp, MTP_results@adjp)
colnames(MTP_summary)=c("gene", "alexa", paste(class1,"mean",sep="_"), paste(class2,"mean",sep="_"), "fold_change", "rawp", "adjp")

#For all 'significant' probes
siggene_summary=MTP_summary[MTP_summary[,"rawp"]<0.05,]


#Write final results to file. Choose the right output file for the comparison in question.
write.table(siggene_summary, file=paste(class1,"vs",class2,"siggenes.txt", sep="_"), sep="\t", row.names=FALSE, quote=FALSE)

###To create a heatmap of all data###
setwd("/home/obig/Projects/Colon_Cancer_Resistance/analysis/Sens_vs_Res")
x=as.matrix(log2(filt_data[,1:22]+16))
#col_names=c("MIP101_Rep1","MIP5FUR_Rep1","MIP101_Rep2","MIP5FUR_Rep2","MIP101_Rep3","MIP5FUR_Rep3","MIP101_Rep6","MIP5FUR_Rep6","RKO_Rep1","RKO5FUR_Rep1","RKO_Rep1","RKO_Rep2","RKO_Rep3","RKO5FUR_Rep1","RKO5FUR_Rep2","RKO5FUR_Rep3","HCT_Rep1","HCT_Rep2","HCT_Rep3","HCT5FUR_Rep1","HCT5FUR_Rep2","HCT5FUR_Rep3")
col_names=c("MIP101_1","MIP5FUR_1","MIP101_2","MIP5FUR_2","MIP101_3","MIP5FUR_3","MIP101_6","MIP5FUR_6","RKO_1","RKO5FUR_1","RKO_1","RKO_2","RKO_3","RKO5FUR_1","RKO5FUR_2","RKO5FUR_3","HCT_1","HCT_2","HCT_3","HCT5FUR_1","HCT5FUR_2","HCT5FUR_3")
cond_colors=c("MIP101","MIP5FUR","MIP101","MIP5FUR","MIP101","MIP5FUR","MIP101","MIP5FUR","RKO","RKO5FUR","RKO","RKO","RKO","RKO5FUR","RKO5FUR","RKO5FUR","HCT","HCT","HCT","HCT5FUR","HCT5FUR","HCT5FUR")
cond_colors[cond_colors=="MIP101"]="#bbbddc"
cond_colors[cond_colors=="MIP5FUR"]="#756bb1"
cond_colors[cond_colors=="RKO"]="#a1d99b"
cond_colors[cond_colors=="RKO5FUR"]="#31a354"
cond_colors[cond_colors=="HCT"]="#fc9272"
cond_colors[cond_colors=="HCT5FUR"]="#de2d26"


pdf("Sensitive_vs_Resistant_heatmap.pdf")
mypalette<-brewer.pal(9,"Blues")
heatmap.2(x, na.rm = TRUE, scale="none", symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=col_names, col=mypalette, cexCol=0.8, ColSideColors=cond_colors)
dev.off()


#Things to consider?
 #Better color schemes.
 #Sidebar colors for cell type or sens/resistant?
 #Only use genes diff expressed in one or more comparison?






#To perform standard MU test
#MU_test_results = array(0, dimnames = list(var_names_filt, c("pvalue")), dim=c(length(var_names_filt),1))
#MU_test_data=exprs(filt_gcrmaData)

#for (i in 1:length(var_names_filt)){
#   mt_values=MU_test_data[i,c(2,12,16,17,18,19,24)]
#   wt_values=MU_test_data[i,c(1,3,4,5,6,7,8,9,10,11,13,14,15,20,21,22,23)]
#   wilcox_result=wilcox.test(x=mt_values, y=wt_values, alternative="two.sided", paired=FALSE)
#   MU_test_results[i,"pvalue"]=wilcox_result$p.value
#}
#Perform simple multiple testing correction
#pvalues=as.numeric(MU_test_results[,"pvalue"])
#pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
#pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]

