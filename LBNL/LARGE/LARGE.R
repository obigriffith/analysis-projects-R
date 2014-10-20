library(multtest)
working_dir="C:/Users/Obi/Documents/My Dropbox/Projects/LARGE/"
#working_dir="C:/Users/Obi/Documents/My Dropbox/Projects/LARGE/laminin_receptor_genes/"
#working_dir="C:/Users/Obi/Documents/My Dropbox/Projects/LARGE/DG_modifying_genes/"

setwd(working_dir)

#Genes of interest in project: LARGE, POMT1, POMT2, Fukutin (aka. FKTN), FKRP, POMGNT1, LARGE2 (aka. GYLTL1B) and ß3GNT1 (aka. B3GNT1)
#choose one
gene="LARGE"
##DG modifying genes##
#gene="POMT1"
#gene="POMT2"
#gene="FKTN"
#gene="FKRP"
#gene="POMGNT1"
#gene="GYLTL1B"
#gene="B3GNT1"
##laminin receptor genes##
#gene="DAG1" #aka DAG, Dystroglycan
#gene="SDC1" 
#gene="RPL10"
#gene="BCAM"
#gene="ITGB1"
#gene="ITGB4"
#gene="ITGA3"
#gene="ITGA1"
#gene="ITGA2"
#gene="ITGA6"
#gene="ITGA7"
#gene="ITGA8"
#gene="ITGA9"
#gene="ITGA5"
#gene="UGT8" # aka Ugt8a
#gene="GAL3ST1"

#Cell line info
celllinefile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/cell_line_info.2.txt"

#All gene expression data
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_cell_data=read.table(celllinefile, header = TRUE, na.strings = "NA", sep="\t")
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Set output files
gene_expression_all_lines_file=paste(gene,"_expression_all_lines.pdf",sep="")
gene_expression_all_lines_file2=paste(gene,"_expression_barplot_all_lines.pdf",sep="")

gene_subtype_comp_file=paste(gene,"_subtype_comp.pdf",sep="")
gene_ERBB2amp_comp_file=paste(gene,"_ERBB2amp_comp.pdf",sep="")
gene_stat_results_file=paste(gene,"_wilcox_test_results.txt",sep="")
Allgenes_DE_results_file="Allgenes_DE_results.txt"

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:length(colnames(raw_data_import))]
raw_exp_status=raw_exp_status_import[,5:length(colnames(raw_data_import))]

#Make sure that cell line info and raw data are ordered the same!!!
libs=colnames(raw_data)
lib_names=as.vector(raw_cell_data[,"Sample.Name"])
cbind(lib_names,libs)

#If ok, use clean names for data object
colnames(raw_data)=lib_names
colnames(raw_exp_status)=lib_names

#Exclude low quality and "Unknown/Normal subtype" libs
#high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"SubtypeNew"]!="Unknown") & as.vector(raw_cell_data[,"SubtypeNew"]!="Normal"))

#Or, include only cell lines also available
#Note: No RNAseq data for MDAMB468, MDAMB436
#Note: MDAMB415 is low quality
#Note: Subtype unknown for HCC1954, MDAMB175VII, SUM159PT in new subtype classification. Use old determinations.
#Note: HCC1500 - identity called into question. But already included throughout publication. Leave for now?

LARGE_study_cell_lines=c("AU565","BT474","CAMA1","LY2","MCF7","MDAMB134VI","MDAMB175VII","MDAMB361","MDAMB415","SKBR3","SUM52PE","T47D","UACC893","ZR75B","HCC70","HCC1143","HCC1419","HCC1569","HCC1937","HCC1954","MDAMB468","BT549","HCC38","HCC1500","HS578T","MDAMB231","MDAMB436","SUM149PT","SUM159PT")
#LARGE_study_cell_lines=c("BT549","HCC70","HCC1143","HCC1500","HCC1569","HCC1937","SUM149PT","AU565","BT474","HCC38","HS578T","MDAMB231","CAMA1","HCC1419","HCC1954","LY2","MCF7","MDAMB134VI","MDAMB175VII","MDAMB361","SKBR3","SUM52PE","SUM159PT","T47D","UACC893","ZR75B")

#high_qual=which(raw_cell_data[,"Sample.Name"] %in% LARGE_study_cell_lines & (as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int"))
high_qual=which(raw_cell_data[,"Sample.Name"] %in% LARGE_study_cell_lines & (as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int" | as.vector(raw_cell_data[,"Sample.Name"])=="HCC1500"))

#Apply library filter to datasets
data=raw_data[,high_qual]
exp_status=raw_exp_status[,high_qual]
cell_data=raw_cell_data[high_qual,]
feat_data=raw_feat_data
libs=colnames(data)

#Retrieve cell line details
lib_names=as.vector(cell_data[,"Sample.Name"])
#subtypes=as.vector(cell_data[,"SubtypeNew"])
subtypes=as.vector(cell_data[,"subtype"]) #Use old subtype definitions for subset of LARGE study lines
ERBB2=as.vector(cell_data[,"ERBB2New"])
Qualities=as.vector(cell_data[,"Quality"])
cbind(lib_names,libs)

#Sort cell lines by subtype then, name
lib_names_sorted=lib_names[order(cell_data[,"subtype"], cell_data[,"Sample.Name"])]
#lib_names_sorted=lib_names[order(cell_data[,"SubtypeNew"], cell_data[,"Sample.Name"])]
subtypes_sorted=subtypes[order(subtypes)]

#Do same thing for HER2+ versus HER2-
#Sort cell lines by ERBB2 then, then name
lib_names_sorted2=lib_names[order(cell_data[,"ERBB2New"], cell_data[,"Sample.Name"])]
ERBB2_sorted=ERBB2[order(ERBB2)]

#Get gene data (LARGE and any other genes of interest), and sort according to subtype order determined above
GENE_data_sorted=data[which(feat_data[,"Seq_Name"]==gene),lib_names_sorted]

#Also get expression status for these genes and summarize
GENE_exp_status_sorted=exp_status[which(feat_data[,"Seq_Name"]==gene),lib_names_sorted]
num_basal=length(which(subtypes_sorted=="Basal"))
num_claudinlow=length(which(subtypes_sorted=="ClaudinLow"))
num_luminal=length(which(subtypes_sorted=="Luminal"))
num_basal_exp=sum(GENE_exp_status_sorted[subtypes_sorted=="Basal"])
num_claudinlow_exp=sum(GENE_exp_status_sorted[subtypes_sorted=="ClaudinLow"])
num_luminal_exp=sum(GENE_exp_status_sorted[subtypes_sorted=="Luminal"])

#Convert values to log2, unless they are 0 in which case set them to 0
GENE_data_sorted_log2 = log2(GENE_data_sorted+1)

#Produce box and whisker plot
subtype_rainbow=rainbow(length(unique(subtypes_sorted)))
names(subtype_rainbow)=unique(subtypes_sorted)
subtype_colors=subtype_rainbow[subtypes_sorted]

#Replace zero values with NA. Works better with boxplot (and is consistent with alexa-seq plots)
data_noZero=data[,lib_names_sorted]
data_noZero[data_noZero==0]=NA
GENE_data_noZero=GENE_data_sorted
GENE_data_noZero[GENE_data_noZero==0]=NA

#Fix cell line names for consistency with LARGE paper
colnames(data_noZero)[colnames(data_noZero)=="MDAMB231"]="MDA231"
colnames(data_noZero)[colnames(data_noZero)=="MDAMB134VI"]="MDA134VI"
colnames(data_noZero)[colnames(data_noZero)=="MDAMB175VII"]="MDA175VII"
colnames(data_noZero)[colnames(data_noZero)=="MDAMB361"]="MDA361"

#Histogram of gene expression data showing position of current gene
pdf(file=gene_expression_all_lines_file, width=10, height=7.5)
par(oma=c(1,1,1,1), mar=c(7,5,2,1)) #bottom, left, top, right
#main_title = "LARGE gene-level expression compared to distribution of all genes"
y_label = "Log2 gene expression level"
boxplot(x = log2(data_noZero), col=subtype_colors, ylab=y_label, las=2, col.lab = gray(.1), cex.main = 2, cex.lab = 1.7, cex.axis=1.5, main=gene)
target_gene_level = as.numeric(log2(GENE_data_noZero))
points(target_gene_level, pch=16, col="black", cex=2)
#legend("bottomleft", legend=c(unique(subtypes_sorted),gene), fill=c(unique(subtype_colors),"white"), pch=c(NA,NA,NA,16), bg="white", border=NA)
dev.off()

#Barplot of gene expression data
barplot_data=as.numeric(GENE_data_sorted)
names(barplot_data)=colnames(GENE_data_sorted)
ylim_value=max(barplot_data)+50
pdf(file=gene_expression_all_lines_file2, width=10, height=7.5)
par(oma=c(1,1,1,1), mar=c(7,5,2,1)) #bottom, left, top, right
barplot(barplot_data, col=subtype_colors, ylab="LARGE Expression", las=2, ylim=c(0,ylim_value), cex.lab = 1.4, cex.axis=1.2, cex.names=1.2)
dev.off()

#Calculate statistics for LARGE differential expression between subtypes and ERBB2 status
###Subtype statistics###
#GENE_subtype_data_log2=list(Basal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Basal"]), Basal_NM=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Basal_NM"]), ClaudinLow=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Luminal"])) 
GENE_subtype_data_log2=list(Basal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Basal"]), ClaudinLow=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Luminal"]))
GENE_subtype_data=list(Basal=as.numeric(GENE_data_sorted[subtypes_sorted=="Basal"]), ClaudinLow=as.numeric(GENE_data_sorted[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(GENE_data_sorted[subtypes_sorted=="Luminal"]))

#Subtype statistics
LvsBC.wilcox=wilcox.test(x=GENE_subtype_data_log2$Luminal, y=c(GENE_subtype_data_log2$Basal,GENE_subtype_data_log2$ClaudinLow), alternative="two.sided")
BvsCL.wilcox=wilcox.test(x=GENE_subtype_data_log2$Basal, y=c(GENE_subtype_data_log2$ClaudinLow,GENE_subtype_data_log2$Luminal), alternative="two.sided")
CvsBL.wilcox=wilcox.test(x=GENE_subtype_data_log2$ClaudinLow, y=c(GENE_subtype_data_log2$Basal,GENE_subtype_data_log2$Luminal), alternative="two.sided")
BvsC.wilcox=wilcox.test(x=GENE_subtype_data_log2$Basal, y=GENE_subtype_data_log2$ClaudinLow, alternative="two.sided")
BvsL.wilcox=wilcox.test(x=GENE_subtype_data_log2$Basal, y=GENE_subtype_data_log2$Luminal, alternative="two.sided")
CvsL.wilcox=wilcox.test(x=GENE_subtype_data_log2$ClaudinLow, y=GENE_subtype_data_log2$Luminal, alternative="two.sided")

#Fold-changes
LvsBC.FC=mean(GENE_subtype_data$Luminal)/mean(c(GENE_subtype_data$Basal,GENE_subtype_data$ClaudinLow))
BvsCL.FC=mean(GENE_subtype_data$Basal)/mean(c(GENE_subtype_data$ClaudinLow,GENE_subtype_data$Luminal))
CvsBL.FC=mean(GENE_subtype_data$ClaudinLow)/mean(c(GENE_subtype_data$Basal,GENE_subtype_data$Luminal))
BvsC.FC=mean(GENE_subtype_data$Basal)/mean(GENE_subtype_data$ClaudinLow)
BvsL.FC=mean(GENE_subtype_data$Basal)/mean(GENE_subtype_data$Luminal)
CvsL.FC=mean(GENE_subtype_data$ClaudinLow)/mean(GENE_subtype_data$Luminal)

pdf(file=gene_subtype_comp_file, width=7.5, height=7.5)
par(oma=c(1,1,1,1), mar=c(7,5,2,1)) #bottom, left, top, right
boxplot(x=GENE_subtype_data_log2, las=2, col=subtype_rainbow, ylab=y_label, cex.axis=1.5, cex.lab=1.5, cex.main=1.7, main=gene)
dev.off()

###ERBB2 statistics (ERBB2pos vs ERBB2neg)
GENE_data_sorted2=data[which(feat_data[,"Seq_Name"]==gene),lib_names_sorted2]
GENE_data_sorted2_log2 = log2(GENE_data_sorted2+1)

GENE_ERBB2_data_log2=list(ERBB2pos=as.numeric(GENE_data_sorted2_log2[ERBB2_sorted=="Amp"]), ERBB2neg=as.numeric(GENE_data_sorted2_log2[ERBB2_sorted=="NoAmp"])) 
GENE_ERBB2_data=list(ERBB2pos=as.numeric(GENE_data_sorted2[ERBB2_sorted=="Amp"]), ERBB2neg=as.numeric(GENE_data_sorted2[ERBB2_sorted=="NoAmp"])) 

HER2posvsHER2neg.wilcox=wilcox.test(x=GENE_ERBB2_data_log2$ERBB2pos, y=GENE_ERBB2_data_log2$ERBB2neg, alternative="two.sided")
HER2posvsHER2neg.FC=mean(GENE_ERBB2_data$ERBB2pos)/mean(GENE_ERBB2_data$ERBB2neg)

pdf(file=gene_ERBB2amp_comp_file, width=7.5, height=7.5)
par(oma=c(1,1,1,1), mar=c(7,5,2,1)) #bottom, left, top, right
boxplot(x=GENE_ERBB2_data_log2, las=2, col=subtype_rainbow, ylab=y_label, cex.axis=1.5, cex.lab=1.5, cex.main=1.7, main=gene)
dev.off()

#Summarize statistics and write to file
tests=c("Luminal vs Basal/ClaudinLow", "Basal vs ClaudinLow/Luminal", "ClaudinLow vs Basal/Luminal", "Basal vs ClaudinLow", "Basal vs Luminal", "ClaudinLow vs Luminal", "HER2+ vs HER2-")
pvalues=c(LvsBC.wilcox$p.value, BvsCL.wilcox$p.value, CvsBL.wilcox$p.value,BvsC.wilcox$p.value,BvsL.wilcox$p.value,CvsL.wilcox$p.value,HER2posvsHER2neg.wilcox$p.value)
FCs=c(LvsBC.FC, BvsCL.FC, CvsBL.FC,BvsC.FC,BvsL.FC,CvsL.FC,HER2posvsHER2neg.FC)

wilcox_results=data.frame(cbind(gene=gene, test=tests, p=pvalues, FC=FCs))
write.table(wilcox_results,file=gene_stat_results_file, sep="\t", row.names=FALSE)

expressed_results=data.frame(cbind(gene=gene, subtype=c("Basal","ClaudinLow","Luminal"), number_expressed=c(num_basal_exp,num_claudinlow_exp,num_luminal_exp), total=c(num_basal,num_claudinlow,num_luminal)))
write.table(expressed_results,file=gene_stat_results_file, sep="\t", row.names=FALSE, append=TRUE)

#Calculate p-value and FC for all genes
###Note, commented out after run once, since for loop is very slow###
#Luminal vs Basal/Claudin-Low
#genes=feat_data[,"EnsEMBL_Gene_ID"]
#gene_names=feat_data[,"Seq_Name"]

#DE_gene_results=matrix(NA,nrow=length(genes),ncol=4, dimnames = list(genes,c("Gene","Symbol","wilcox.pvalue","FC")))
#for (i in 1:length(genes)){
# gene=genes[i]
# gene_name=gene_names[i]
# GENE_data_sorted=data[which(feat_data[,"EnsEMBL_Gene_ID"]==gene),lib_names_sorted]
# GENE_data_sorted_log2 = log2(GENE_data_sorted+1)
# GENE_subtype_data_log2=list(Basal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Basal"]), ClaudinLow=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(GENE_data_sorted_log2[subtypes_sorted=="Luminal"]))
# GENE_subtype_data=list(Basal=as.numeric(GENE_data_sorted[subtypes_sorted=="Basal"]), ClaudinLow=as.numeric(GENE_data_sorted[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(GENE_data_sorted[subtypes_sorted=="Luminal"]))
# LvsBC.wilcox=wilcox.test(x=GENE_subtype_data_log2$Luminal, y=c(GENE_subtype_data_log2$Basal,GENE_subtype_data_log2$ClaudinLow), alternative="two.sided")
# LvsBC.FC=mean(GENE_subtype_data$Luminal)/mean(c(GENE_subtype_data$Basal,GENE_subtype_data$ClaudinLow))
# DE_gene_results[i,"Gene"]=gene
# DE_gene_results[i,"Symbol"]=gene_name
# DE_gene_results[i,"wilcox.pvalue"]=LvsBC.wilcox$p.value
# DE_gene_results[i,"FC"]=LvsBC.FC
#}

#Correct p-values
#pvalues=as.numeric(DE_gene_results[,"wilcox.pvalue"])
#pvalues_adj=mt.rawp2adjp(pvalues, proc=c("BH"))
#pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
#DE_gene_results=cbind(DE_gene_results, pvalues_adj_orig_order[,2])

#write.table(DE_gene_results,file=Allgenes_DE_results_file, sep="\t", row.names=FALSE)
