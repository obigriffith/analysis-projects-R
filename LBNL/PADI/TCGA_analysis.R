library("ggplot2")

#specify data dir/files
datadir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/"
setwd(datadir)
TCGA_agilent_normal_datafile="AgilentExpression_V2/TCGA_BrCa_AgilentG4502A_NormalsMatched_March72012_KNNImputed2_V2.txt"
TCGA_agilent_tumor_datafile="AgilentExpression_V2/TCGA_BrCa_AgilentG4502A_Tumors_KNNImputed2_V2.txt"
TCGA_rnaseq_normal_datafile="RNAseq_V2/TCGA_BrCa_IlluminaHiSeq_RNASeq_NormalsMatched_March72012_Postprocessed_Gene_RPKM.csv"
TCGA_rnaseq_tumor_datafile="RNAseq_V2/TCGA_BrCa_IlluminaHiSeq_RNASeq_Tumors_March72012_Postprocessed_Gene_RPKM.csv"

#Specify gene of interest
gene="PADI4"
RSgene="PADI4|23569" #Absurd gene ID from TCGA RNAseq data
RSgene2="GSK3B|2932"

#Load data
TCGA_agilent_normal_import=read.table(TCGA_agilent_normal_datafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)
TCGA_agilent_tumor_import=read.table(TCGA_agilent_tumor_datafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)
TCGA_rnaseq_normal_import=read.csv(TCGA_rnaseq_normal_datafile, header = TRUE, na.strings = "NA", row.names=1)
TCGA_rnaseq_tumor_import=read.csv(TCGA_rnaseq_tumor_datafile, header = TRUE, na.strings = "NA", row.names=1)

#Set result dir
resultdir="C:/Users/Obi/Documents/My Dropbox/Projects/PADI/TCGA_analysis/"
setwd(resultdir)

#Make sample names consistent to allow merging of tumor/normal data
colnames(TCGA_agilent_normal_import)=substr(colnames(TCGA_agilent_normal_import),1,12)
colnames(TCGA_agilent_tumor_import)=substr(colnames(TCGA_agilent_tumor_import),1,12)
agilent_normal_samples=colnames(TCGA_agilent_normal_import)
agilent_tumor_samples=colnames(TCGA_agilent_tumor_import)
agilent_matched_samples=agilent_normal_samples[which(agilent_normal_samples %in% agilent_tumor_samples)]

colnames(TCGA_rnaseq_normal_import)=substr(colnames(TCGA_rnaseq_normal_import),1,12)
colnames(TCGA_rnaseq_tumor_import)=substr(colnames(TCGA_rnaseq_tumor_import),1,12)
rnaseq_normal_samples=colnames(TCGA_rnaseq_normal_import)
rnaseq_tumor_samples=colnames(TCGA_rnaseq_tumor_import)
rnaseq_matched_samples=rnaseq_normal_samples[which(rnaseq_normal_samples %in% rnaseq_tumor_samples)]
rnaseq_matched_samples=rnaseq_matched_samples[-which(rnaseq_matched_samples=="X")] #Exclude "X" column

#Extract data for gene and patients of interest
agilent_normal_gene_data=TCGA_agilent_normal_import[gene,agilent_matched_samples]
agilent_tumor_gene_data=TCGA_agilent_tumor_import[gene,agilent_matched_samples]

#PADI4
rnaseq_normal_gene_data=TCGA_rnaseq_normal_import[RSgene,rnaseq_matched_samples]
rnaseq_tumor_gene_data=TCGA_rnaseq_tumor_import[RSgene,rnaseq_matched_samples]

#GSK3B
rnaseq_normal_gene2_data=TCGA_rnaseq_normal_import[RSgene2,rnaseq_matched_samples]
rnaseq_tumor_gene2_data=TCGA_rnaseq_tumor_import[RSgene2,rnaseq_matched_samples]

#Test for difference between normal/tumor for gene of interest
agilent.wilcox.pvalue=wilcox.test(x=as.numeric(agilent_normal_gene_data), y=as.numeric(agilent_tumor_gene_data), alternative="two.sided", paired=TRUE)$p.value
rnaseq.wilcox.pvalue=wilcox.test(x=as.numeric(rnaseq_normal_gene_data), y=as.numeric(rnaseq_tumor_gene_data), alternative="two.sided", paired=TRUE)$p.value
rnaseq.wilcox.pvalue2=wilcox.test(x=as.numeric(rnaseq_normal_gene2_data), y=as.numeric(rnaseq_tumor_gene2_data), alternative="two.sided", paired=TRUE)$p.value

#Create dataframe appropriate for ggplots
ExpVal=c(as.numeric(agilent_normal_gene_data),as.numeric(agilent_tumor_gene_data))
Type=c(rep("normal",length(agilent_matched_samples)),rep("tumor",length(agilent_matched_samples)))
agilent_gene_data=data.frame(ExpVal=ExpVal,Type=Type)

ExpVal=c(as.numeric(rnaseq_normal_gene_data),as.numeric(rnaseq_tumor_gene_data))
Type=c(rep("normal",length(rnaseq_matched_samples)),rep("tumor",length(rnaseq_matched_samples)))
rnaseq_gene_data=data.frame(ExpVal=ExpVal,Type=Type)

ExpVal=c(as.numeric(rnaseq_normal_gene2_data),as.numeric(rnaseq_tumor_gene2_data))
Type=c(rep("normal",length(rnaseq_matched_samples)),rep("tumor",length(rnaseq_matched_samples)))
rnaseq_gene2_data=data.frame(ExpVal=ExpVal,Type=Type)


#Create plots showing distributions
pdf(file="TCGA_Agilent_PADI4.pdf")
ggplot(agilent_gene_data, aes(ExpVal, fill=Type)) + geom_density(alpha = 0.2, stat="bin", binwidth=( max(agilent_gene_data[,"ExpVal"]) - min(agilent_gene_data[,"ExpVal"]) )/30 ) + xlab("Expression value") + ylab("Count") + opts(title=gene)
dev.off()

pdf(file="TCGA_RNAseq_PADI4.pdf")
ggplot(rnaseq_gene_data, aes(ExpVal, fill=Type)) + geom_density(alpha = 0.2, stat="bin", binwidth=( max(rnaseq_gene_data[,"ExpVal"]) - min(rnaseq_gene_data[,"ExpVal"]) )/30 ) + xlab("Expression value") + ylab("Count") + opts(title=RSgene)
dev.off()

pdf(file="TCGA_RNAseq_GSK3B.pdf")
ggplot(rnaseq_gene2_data, aes(ExpVal, fill=Type)) + geom_density(alpha = 0.2, stat="bin", binwidth=( max(rnaseq_gene2_data[,"ExpVal"]) - min(rnaseq_gene2_data[,"ExpVal"]) )/30 ) + xlab("Expression value") + ylab("Count") + opts(title=RSgene2)
dev.off()


#Agilent - Boxplot of PADI4 by type (tumor vs normal)
agilent_gene_subtype_data=list(Normal=as.numeric(agilent_normal_gene_data), Tumor=as.numeric(agilent_tumor_gene_data)) 
pdf(file="PADI4_Agilent_normal_vs_tumor_comp.pdf", width=7.5, height=7.5)
y_label = "Agilent Expression Value"
subtype_rainbow=rainbow(2)
ymin=min(as.numeric(agilent_gene_data[,"ExpVal"]))
ymax=max(as.numeric(agilent_gene_data[,"ExpVal"]))
#ymax=1
par(mar=c(8,5,3,1)) #bottom, left, top, right
boxplot(x=agilent_gene_subtype_data, las=2, col=subtype_rainbow, main="PADI4", ylab=y_label, ylim=c(ymin,ymax), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
legend("topright",legend=paste("p = ",format(agilent.wilcox.pvalue, digits=3),sep=""),bty="n")
dev.off()

#RNAseq - Boxplot of PADI4 by type (tumor vs normal)
rnaseq_gene_subtype_data=list(Normal=as.numeric(log2(rnaseq_normal_gene_data+1)), Tumor=as.numeric(log2(rnaseq_tumor_gene_data+1))) 
pdf(file="PADI4_RNAseq_normal_vs_tumor_comp.pdf", width=7.5, height=7.5)
y_label = "RNA-seq log2(FPKM+1)"
subtype_rainbow=rainbow(2)
ymin=min(as.numeric(rnaseq_gene_data[,"ExpVal"]))
#ymax=max(as.numeric(rnaseq_gene_data[,"ExpVal"]))
ymax=0.8
par(mar=c(8,5,3,1)) #bottom, left, top, right
boxplot(x=rnaseq_gene_subtype_data, las=2, col=subtype_rainbow, main="PADI4", ylab=y_label, ylim=c(ymin,ymax), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
legend("topright",legend=paste("p = ",format(rnaseq.wilcox.pvalue, digits=3),sep=""),bty="n")
dev.off()


#RNAseq - Boxplot of GSK3B by type (tumor vs normal)
rnaseq_gene2_subtype_data=list(Normal=as.numeric(log2(rnaseq_normal_gene2_data+1)), Tumor=as.numeric(log2(rnaseq_tumor_gene2_data+1))) 
pdf(file="GSK3B_RNAseq_normal_vs_tumor_comp.pdf", width=7.5, height=7.5)
y_label = "RNA-seq log2(FPKM+1)"
subtype_rainbow=rainbow(2)
#ymin=min(as.numeric(rnaseq_gene2_data[,"ExpVal"]))
ymin=0
#ymax=max(as.numeric(rnaseq_gene2_data[,"ExpVal"]))
ymax=5
par(mar=c(8,5,3,1)) #bottom, left, top, right
boxplot(x=rnaseq_gene2_subtype_data, las=2, col=subtype_rainbow, main="GSK3B", ylab=y_label, ylim=c(ymin,ymax), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
legend("topright",legend=paste("p = ",format(rnaseq.wilcox.pvalue2, digits=3),sep=""),bty="n")
dev.off()

#Determine correlation between PADI4 and GSK3B in RNAseq data
PADI4_GSK3B_rnaseq_normal_rho=cor.test(x=as.numeric(rnaseq_normal_gene_data), y=as.numeric(rnaseq_normal_gene2_data), method ="spearman")
PADI4_GSK3B_rnaseq_tumor_rho=cor.test(x=as.numeric(rnaseq_tumor_gene_data), y=as.numeric(rnaseq_tumor_gene2_data), method ="spearman")
PADI4_GSK3B_rnaseq_all_rho=cor.test(x=c(as.numeric(rnaseq_normal_gene_data),as.numeric(rnaseq_tumor_gene_data)), y=c(as.numeric(rnaseq_normal_gene2_data),as.numeric(rnaseq_tumor_gene2_data)), method ="spearman")



pdf(file="PADI4_GSK3B_cor_by_type.pdf", width=10, height=7.5)
#Break plot into sections for each subtype to be plotted proportionally
layout(matrix(c(1,1,1,1,1,2,2,2,2,2), 1, 10, byrow = TRUE))
#layout.show(2)

par(mar=c(7,4,3,1)) #bottom, left, top, right
#Create plot to represent correlation between PADI2 and GSK3B
#break into panels for normal/tumor
PADI4_data_log2_normal=log2(as.numeric(rnaseq_normal_gene_data)+1)
GSK3B_data_log2_normal=log2(as.numeric(rnaseq_normal_gene2_data)+1)
plot(x=1:length(PADI4_data_log2_normal), y=PADI4_data_log2_normal, type="n", xlim=c(1,length(PADI4_data_log2_normal)), ylim=c(0,max(c(PADI4_data_log2_normal,GSK3B_data_log2_normal,PADI4_data_log2_tumor,GSK3B_data_log2_tumor))), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", bty="n")
axis(1, at=1:length(PADI4_data_log2_normal), labels=NA, las=2, cex.axis=0.9)
points(x=(1:length(PADI4_data_log2_normal)), y=PADI4_data_log2_normal, col="blue", pch=18, cex=1.6)
points(x=(1:length(GSK3B_data_log2_normal)), y=GSK3B_data_log2_normal, col="red", pch=15, cex=1.2)
lines(x=(1:length(PADI4_data_log2_normal)), y=PADI4_data_log2_normal, col="blue", lwd=2)
lines(x=(1:length(GSK3B_data_log2_normal)), y=GSK3B_data_log2_normal, col="red", lwd=2)
legend("top", legend=c("Normal", paste("rho =",format(PADI4_GSK3B_rnaseq_normal_rho$estimate, digit=3)),paste("p =",format(PADI4_GSK3B_rnaseq_normal_rho$p.value,digit=3))), bty="n")
legend("topleft", legend=c("PADI4","GSK3B"), pch=c(18,15), lwd=c(2,2), col=c("blue","red"), pt.cex=c(1.5,1.2), bty="n")

par(mar=c(7,2,3,1)) #bottom, left, top, right
PADI4_data_log2_tumor=log2(as.numeric(rnaseq_tumor_gene_data)+1)
GSK3B_data_log2_tumor=log2(as.numeric(rnaseq_tumor_gene2_data)+1)
plot(x=1:length(PADI4_data_log2_tumor), y=PADI4_data_log2_tumor, type="n", xlim=c(1,length(PADI4_data_log2_tumor)), ylim=c(0,max(c(PADI4_data_log2_normal,GSK3B_data_log2_normal,PADI4_data_log2_tumor,GSK3B_data_log2_tumor))), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", yaxt="n", bty="n")
axis(1, at=1:length(PADI4_data_log2_tumor), labels=NA, las=2, cex.axis=0.9)
points(x=(1:length(PADI4_data_log2_tumor)), y=PADI4_data_log2_tumor, col="blue", pch=18, cex=1.6)
points(x=(1:length(GSK3B_data_log2_tumor)), y=GSK3B_data_log2_tumor, col="red", pch=15, cex=1.2)
lines(x=(1:length(PADI4_data_log2_tumor)), y=PADI4_data_log2_tumor, col="blue", lwd=2)
lines(x=(1:length(GSK3B_data_log2_tumor)), y=GSK3B_data_log2_tumor, col="red", lwd=2)
legend("top", legend=c("Tumor", paste("rho =",format(PADI4_GSK3B_rnaseq_tumor_rho$estimate, digit=3)),paste("p =",format(PADI4_GSK3B_rnaseq_tumor_rho$p.value,digit=3))), bty="n")
dev.off()


