#Load the appropriate libraries
library(affy)
library(gcrma)
library(multtest)
library(Biobase)
library(genefilter)
library("biomaRt")

#Set working directory for data files
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/CEL_FILES/tumor_GC")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/CEL_FILES/tumor")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/CEL_FILES/tumor_plus_unsorted")
setwd("/home/obig/Projects/Lymphoma_drugs/FL/CEL_FILES/tumor_EZH2_status")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/CEL_FILES/tumor_micro")
#setwd("/home/obig/Projects/Lymphoma_drugs/DLBCL/CEL_FILES/")

Data=ReadAffy()

#Run GCRMA on Data to background correct, normalize and summarize expression data
gcrmaData=gcrma(Data)

#Set working data for results files
setwd("/home/obig/Projects/Lymphoma_drugs/FL")
#setwd("/home/obig/Projects/Lymphoma_drugs/DLBCL")

#Write gcrma normalized data to file
#First, reduce number of decimal places
gcrma_values_formatted=format(exprs(gcrmaData), digits=5)
#write.table(gcrma_values_formatted, file = "FL_vs_normGCBcell_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "FL_vs_normGCBcell_gcrma.2.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "DLBCL_vs_normGCBcell_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "FL_vs_micro_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "FL_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
#write.table(gcrma_values_formatted, file = "FL_gcrma.2.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
write.table(gcrma_values_formatted, file = "FL_EZH2_status_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")

#Set working data for results files
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/up_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/down_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/diff_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/FL/tumor_vs_micro/diff_regulated")
setwd("/home/obig/Projects/Lymphoma_drugs/FL/EZH2_mt_vs_wt/diff_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/DLBCL/up_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/DLBCL/down_regulated")
#setwd("/home/obig/Projects/Lymphoma_drugs/DLBCL/diff_regulated")

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

#Create a vector of probe names after filtering for output later
var_names_filt=rownames(exprs(filt_gcrmaData))

#Then perform differential expression statistics and multiple testing correction
#Create a vector for the class labels ('FL'=0 or 'DLBCL'=0 vs 'GC B-cell'=1).  This determines what values are tested against what

#FL Data:
#class_labels=vector(length=29)
#class_labels[1:24]="Tumor"
#class_labels[25:29]="Normal"

#tumor vs micro
#class_labels=vector(length=48)
#class_labels[c(2,4,6,8,10,12,14,16,18,20,22,24,26,27,29,31,33,35,37,39,41,43,45,48)]="Tumor"
#class_labels[c(1,3,5,7,9,11,13,15,17,19,21,23,25,28,30,32,34,36,38,40,42,44,46,47)]="Micro"

#EZH2 mt vs wt
class_labels=vector(length=24)
class_labels[c(2,12,16,17,18,19,24)]="mutant"
class_labels[c(1,3,4,5,6,7,8,9,10,11,13,14,15,20,21,22,23)]="wildtype"

#DLBCL Data:
#class_labels=vector(length=26)
#class_labels[1:5]="Normal"
#class_labels[6:26]="Tumor"

#Filter data down to just probes with a PCG target
PCG_data_temp=exprs(filt_gcrmaData)
#For some reason all probes can't be specified in a single command. Too many items for a single command? Break into several
PCG_probe_list_A=c("1552626_a_at","1552733_at","1553530_a_at","1553678_a_at","1554690_a_at","1557905_s_at","200632_s_at","200661_at","200764_s_at","200765_x_at","200953_s_at","201236_s_at","201464_x_at","201465_s_at","201466_s_at","202016_at","202156_s_at","202241_at","202284_s_at","202524_s_at","202638_s_at","202726_at","202897_at","202935_s_at","202936_s_at","203037_s_at","203178_at","203232_s_at","203395_s_at","203408_s_at","203574_at","203603_s_at","203684_s_at","203685_at","203688_at","203695_s_at","203773_x_at","203881_s_at","203921_at","204034_at","204285_s_at","204286_s_at","204301_at","204451_at","204480_s_at","204489_s_at","204490_s_at","204491_at","204591_at","204745_x_at","204805_s_at","204897_at","204995_at","205081_at","205128_x_at","205249_at","205466_s_at","205590_at","205632_s_at","205691_at","205758_at","205780_at","205885_s_at","205990_s_at","206045_s_at","206115_at","206194_at","206693_at","206864_s_at","207001_x_at","207087_x_at","207425_s_at","207480_s_at","207522_s_at","207813_s_at")
PCG_probe_list_B=c("207828_s_at","208018_s_at","208322_s_at","208352_x_at","208353_x_at","208657_s_at","208763_s_at","208891_at","208892_s_at","208893_s_at","208999_at","209030_s_at","209135_at","209277_at","209278_s_at","209447_at","209604_s_at","209644_x_at","209681_at","209782_s_at","209824_s_at","209835_x_at","210154_at","210786_s_at","210844_x_at","210916_s_at","211026_s_at","211323_s_at","211729_x_at","211752_s_at","211825_s_at","212014_x_at","212099_at","212188_at","212192_at","212390_at","212509_s_at","212707_s_at","212720_at","212977_at","213036_x_at","213082_s_at","213083_at","213093_at","213147_at","213150_at","213294_at","213353_at","213385_at","213425_at","213497_at","213523_at","213844_at","214157_at","214639_s_at","214933_at","215813_s_at","215990_s_at","216733_s_at","217762_s_at","217763_s_at","217764_s_at","218532_s_at","218600_at","218739_at","218764_at","218866_s_at","218974_at","219359_at","219634_at","221211_s_at","221234_s_at","221584_s_at","221666_s_at","221680_s_at","221841_s_at")
PCG_probe_list_C=c("222165_x_at","222816_s_at","223282_at","223287_s_at","223454_at","223503_at","224847_at","224848_at","224851_at","225102_at","225687_at","225803_at","225868_at","225940_at","225941_at","226010_at","226065_at","226333_at","226425_at","226560_at","227173_s_at","227178_at","227361_at","227641_at","227703_s_at","227935_s_at","228868_x_at","228904_at","229256_at","229553_at","229584_at","229642_at","229844_at","230100_x_at","230425_at","230803_s_at","230884_s_at","231776_at","232080_at","235121_at","235228_at","235521_at","235645_at","235753_at","238605_at","238669_at","243932_at","244623_at","41047_at","44783_s_at")
PCG_probe_list=c(PCG_probe_list_A,PCG_probe_list_B,PCG_probe_list_C)
PCG_data=PCG_data_temp[PCG_probe_list,]
var_names_filt=rownames(PCG_data) 


#The Mann-Whitney test (or Wilcoxen Rank Test)="t.twosamp.equalvar".
#You might also consider using the Welch T-test (t.twosamp.unequalvar).  Are variances equal?
#Also, note that B=1000 (the number of permutations) or higher is recommended.  But, for testing 100 is convenient.
#MTP_results=MTP(X=exprs(filt_gcrmaData), Y=class_labels, na.rm=TRUE, alternative="greater", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP")
#MTP_results=MTP(X=exprs(filt_gcrmaData), Y=class_labels, na.rm=TRUE, alternative="less", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP")
MTP_results=MTP(X=exprs(filt_gcrmaData), Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP")
MTP_results=MTP(X=PCG_data, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP")

#Output results.
MTP_summary=cbind(var_names_filt,MTP_results@statistic,MTP_results@estimate,MTP_results@rawp, MTP_results@adjp, MTP_results@reject)
colnames(MTP_summary)=c("probe", "statistic", "estimate", "rawp", "adjp", "reject")

#For only probes significant after multiple testing correction
sigprobe_summary=MTP_summary[MTP_summary[,"reject"]=="TRUE",]

#For all 'significant' probes
sigprobe_summary=MTP_summary[MTP_summary[,"rawp"]<0.05,]
		

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


#write.table(sigprobe_summary, file="FL_vs_normGCBcell.sigprobes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(sigprobe_summary, file="FL_vs_normGCBcell.sigprobes.adjpvalues.2.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(MTP_summary, file="FL_vs_normGCBcell.allprobes.adjpvalues.2.txt",sep="\t",row.names=FALSE,quote=FALSE)
#write.table(sigprobe_summary, file="FL_vs_micro.sigprobes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(MTP_summary, file="FL_vs_micro.allprobes.adjpvalues.txt",sep="\t",row.names=FALSE,quote=FALSE)
#write.table(sigprobe_summary, file="EZH2_mt_vs_wt.sigprobes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(sigprobe_summary, file="EZH2_mt_vs_wt.sigprobes.adjpvalues.2.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(sigprobe_summary, file="EZH2_mt_vs_wt.sigprobes.adjpvalues.3.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(sigprobe_summary, file="EZH2_mt_vs_wt.sigprobes.rawpvalues.3.txt", sep="\t", row.names=FALSE, quote=FALSE)
#write.table(sigprobe_summary, file="DLBCL_vs_normGCBcell.sigprobes.adjpvalues.txt", sep="\t", row.names=FALSE, quote=FALSE)

####################Annotate Probes to Genes or Proteins################################
probe_ids=sigprobe_summary[,1]
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
#write.table(annotations_notnull_unique_unambig, file = "FL_vs_normGCBcell.sigprobes_2_Biomart_uniprot_swissprot_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "DLBCL_vs_normGCBcell.sigprobes_2_Biomart_uniprot_swissprot_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "FL_vs_normGCBcell.sigprobes_2_Biomart_HGNC_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "FL_vs_micro.sigprobes_2_Biomart_HGNC_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "EZH2_mt_vs_wt.sigprobes_2_Biomart_HGNC_unambig_mappings.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "EZH2_mt_vs_wt.sigprobes_2_Biomart_HGNC_unambig_mappings.2.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
#write.table(annotations_notnull_unique_unambig, file = "EZH2_mt_vs_wt.sigprobes_2_Biomart_HGNC_unambig_mappings.3.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
write.table(annotations_notnull_unique_unambig, file = "EZH2_mt_vs_wt.sigprobes_rawp_2_Biomart_HGNC_unambig_mappings.3.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")

write.table(annotations_notnull_unique_unambig, file = "EZH2_mt_vs_wt.allprobes_2_Biomart_HGNC_unambig_mappings.3.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double")
