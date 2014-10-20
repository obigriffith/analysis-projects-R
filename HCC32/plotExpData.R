library(ggplot2)
library(reshape2)
library(grid)
library(genefilter)
library(gplots)
library(multtest)

#This script compares RMA expression estimates from the Malouf et al paper for a set of 39 samples that includes
#pure FL, mixed FL, standard HCC, and normals
#U133Plus2 probe ids were mapped to ENSG and Ensembl Gene Names using Ensembl BioMart

#To do: Reprocess from raw using GCRMA?
#To do: Look up probe sequences for PRKACA probes to see what part of the gene/transcripts they actually are designed for

setwd("/Users/ogriffit/Dropbox/WashU/Projects/HCC32/Analysis/E-MTAB-1503/")
mappingfile="Biomart_U133Plus2_GeneName_mapping.txt"
#datafile="DNAJB1_PRKACA_expression.txt"
datafile="CIT_CHEF_EXP_RMA_DATA.txt"
clinicaldatafile="E-MTAB-1503.sdrf.txt"

#Read in raw data
rawdata=read.table(file=datafile, header=TRUE, sep="\t", row.names=1)
clindata=read.table(file=clinicaldatafile, header=TRUE, sep="\t", row.names=1, as.is=c(1:31))
mappingdata=read.table(file=mappingfile, header=TRUE, sep="\t", na.strings="", as.is=c(1:6))

#Filter to genes that actually have a probe sets specified
mappingdata=mappingdata[which(!is.na(mappingdata[,"Affy.HG.U133.PLUS.2.probeset"])),]

#Filter to only protein_coding and KNOWN genes
mappingdata=mappingdata[which(mappingdata[,"Gene.Biotype"]=="protein_coding" & mappingdata[,"Status..gene."]=="KNOWN"),]

#Filter to just probesets that map to exactly one gene
redundantprobes=names(which(table(mappingdata[,"Affy.HG.U133.PLUS.2.probeset"])>1))
uniqueprobes=names(which(table(mappingdata[,"Affy.HG.U133.PLUS.2.probeset"])==1))
#test=mappingdata[which(mappingdata[,"Affy.HG.U133.PLUS.2.probeset"]%in%redundantprobes),]
#test[order(test[,"Affy.HG.U133.PLUS.2.probeset"]),]
mappingdata=mappingdata[which(mappingdata[,"Affy.HG.U133.PLUS.2.probeset"]%in%uniqueprobes),]

#since probeset names are now unique they can be used for rownames for convenience
rownames(mappingdata)=mappingdata[,"Affy.HG.U133.PLUS.2.probeset"]

#NOTE: The order of operations matters here:
#For example, some probesets will map to both a protein_coding gene and its homologous pseudogene
#If we want to keep that gene we should filter out non-protein-coding first then check to see if it still maps to more than one gene or not
#If we think that such homologous mappings obfuscate what is actually being measured by the probeset then we should filter for unique mapping first
#Other strategies:
#1. proceed with probesets, see what is associated and then sort out what we think is really being measured later
#2. Switch to customCDF annotations where a concerted effort is made to remap probes to a new set of non-redundant/ambiguous probesets

samples=rownames(clindata)
normnames=clindata[samples,"Normalization.Name"]
diseasetypes=clindata[samples,"Characteristics.disease.state."]

#Rename diseasetypes to shorter name for plotting
diseasetypes[which(diseasetypes=="pure fibrolamellar carcinoma")]="pureFL"
diseasetypes[which(diseasetypes=="mixed fibrolamellar carcinoma")]="mixed"
diseasetypes[which(diseasetypes=="non-cirrhotic Hepatocellular carcinoma")]="HCC"
diseasetypes[which(diseasetypes=="non-tumoral adjacent livers")]="Normal"

#Perform PRKACA-specific analysis
PRKACA_mapping=mappingdata[which(mappingdata[,"Associated.Gene.Name"]=="PRKACA"),]
data=data.frame(normnames,t(rawdata[PRKACA_mapping[,"Affy.HG.U133.PLUS.2.probeset"],normnames]),diseasetypes)
colnames(data)=c("sample","PRKACA_216234_s_at","PRKACA_202801_at","disease_state")

#Set levels of factor to desired display order for ggplot
data$disease_state=factor(data$disease_state, levels=c("pureFL","mixed","HCC","Normal"))

#Set up pdf file for output
pdf(file="PRKACA_RMA_vs_disease_type.pdf")
#Set up grid for multi-panel figure
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#Create boxplots
p1 = ggplot(data, aes(factor(disease_state), PRKACA_216234_s_at)) + geom_boxplot(aes(fill = factor(disease_state))) + ylab("RMA value") + xlab("") + theme(legend.position="none") + labs(title="PRKACA (216234_s_at)")
p2 = ggplot(data, aes(factor(disease_state), PRKACA_202801_at)) + geom_boxplot(aes(fill = factor(disease_state))) + ylab("RMA value") + xlab("") + theme(legend.position="none") + labs(title="PRKACA (202801_at)")
#Arrange plots on grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
dev.off()

#Calculate pairwise stats
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_216234_s_at"],y=data[which(data[,"disease_state"]=="mixed"),"PRKACA_216234_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_216234_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"PRKACA_216234_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_216234_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"PRKACA_216234_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"PRKACA_216234_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"PRKACA_216234_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"PRKACA_216234_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"PRKACA_216234_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_202801_at"],y=data[which(data[,"disease_state"]=="mixed"),"PRKACA_202801_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_202801_at"],y=data[which(data[,"disease_state"]=="HCC"),"PRKACA_202801_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"PRKACA_202801_at"],y=data[which(data[,"disease_state"]=="Normal"),"PRKACA_202801_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"PRKACA_202801_at"],y=data[which(data[,"disease_state"]=="HCC"),"PRKACA_202801_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"PRKACA_202801_at"],y=data[which(data[,"disease_state"]=="Normal"),"PRKACA_202801_at"])


#Perform DNAJB1-specific analysis
DNAJB1_mapping=mappingdata[which(mappingdata[,"Associated.Gene.Name"]=="DNAJB1"),]
data=data.frame(normnames,t(rawdata[DNAJB1_mapping[,"Affy.HG.U133.PLUS.2.probeset"],normnames]),diseasetypes)
colnames(data)=c("sample","DNAJB1_200664_s_at","DNAJB1_200666_s_at","DNAJB1_231556_at","disease_state")

#Set levels of factor to desired display order for ggplot
data$disease_state=factor(data$disease_state, levels=c("pureFL","mixed","HCC","Normal"))

#Set up pdf file for output
pdf(file="DNAJB1_RMA_vs_disease_type.pdf", width=10, height=8)
#Set up grid for multi-panel figure
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#Create boxplots
p1 = ggplot(data, aes(factor(disease_state), DNAJB1_200664_s_at)) + geom_boxplot(aes(fill = factor(disease_state))) + ylab("RMA value") + xlab("") + theme(legend.position="none") + labs(title="DNAJB1 (200664_s_at)")
p2 = ggplot(data, aes(factor(disease_state), DNAJB1_200666_s_at)) + geom_boxplot(aes(fill = factor(disease_state))) + ylab("RMA value") + xlab("") + theme(legend.position="none") + labs(title="DNAJB1 (200666_s_at)")
p3 = ggplot(data, aes(factor(disease_state), DNAJB1_231556_at)) + geom_boxplot(aes(fill = factor(disease_state))) + ylab("RMA value") + xlab("") + theme(legend.position="none") + labs(title="DNAJB1 (231556_at)")

#Arrange plots on grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
print(p3, vp = vplayout(1, 3))
dev.off()

#Calculate pairwise stats
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200664_s_at"],y=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200664_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200664_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_200664_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200664_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_200664_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200664_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_200664_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200664_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_200664_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200666_s_at"],y=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200666_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200666_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_200666_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_200666_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_200666_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200666_s_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_200666_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_200666_s_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_200666_s_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_231556_at"],y=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_231556_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_231556_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_231556_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="pureFL"),"DNAJB1_231556_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_231556_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_231556_at"],y=data[which(data[,"disease_state"]=="HCC"),"DNAJB1_231556_at"])
wilcox.test(x=data[which(data[,"disease_state"]=="mixed"),"DNAJB1_231556_at"],y=data[which(data[,"disease_state"]=="Normal"),"DNAJB1_231556_at"])

#Next we will perform global gene expression analysis
#Reduce expression data to just unique filtered probes from above
expdata=rawdata[rownames(mappingdata),normnames]

PRKACA_DNAJB1_data=expdata[c("216234_s_at","202801_at","200664_s_at","200666_s_at","231556_at"),]

#Filter out genes not expressed in many samples or with too low variance to be useful
ffun=filterfun(pOverA(p = 0.1, A = 100), cv(a = 0.5, b = 10))
filt=genefilter(2^expdata,ffun)
expdatafilt=expdata[filt,]

#Create heatmap to represent global gene expression pattern
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Set up colors for sidebar
colors=rainbow(4)
diseasecolors=diseasetypes
diseasecolors[diseasecolors=="pureFL"]=colors[1]
diseasecolors[diseasecolors=="mixed"]=colors[2]
diseasecolors[diseasecolors=="HCC"]=colors[3]
diseasecolors[diseasecolors=="Normal"]=colors[4]

pdf(file="globalGeneExpHeatmap.pdf")
heatmap.2(as.matrix(expdatafilt), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(4,10), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main="", labCol=samples, labRow=FALSE, ColSideColors=diseasecolors, cexCol=0.8, col=rev(heat.colors(75)))
legend("topright",legend=c("pureFL","mixed","HCC","Normal"),fill=c(colors[1],colors[2],colors[3],colors[4]), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#Perform differential expression analysis to see where PRKACA and DNAJB1 stack up overall and to create a list of DE genes for plotting and pathway analysis
#Compare pureFL to HCC/Normal (exclude mixed because they are truly intermediate)

#COMPARISON #1 pureFL vs HCC/Normal
names(diseasetypes)=samples
names(normnames)=samples
samples_sub=samples[which(diseasetypes%in%c("pureFL","HCC","Normal"))]
normnames_sub=normnames[samples_sub]
diseasetypes_sub=diseasetypes[samples_sub]
expdatafilt_sub=expdatafilt[,normnames_sub]
class_labels=diseasetypes_sub
class_labels[which(class_labels=="pureFL")]=0
class_labels[which(class_labels=="HCC" | class_labels=="Normal")]=1

#Then perform differential expression statistics and multiple testing correction for each comparison
#Note to get non-parametric equivalent of t-test, use robust=TRUE
MTP_results=MTP(X=expdatafilt_sub, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=10000, method="sd.minP", robust=TRUE)

#Calculate basic summary stats (means, fold change, etc)
expdatafilt_sub_nonlog2=2^expdatafilt_sub
class1_means=apply(expdatafilt_sub_nonlog2[,which(class_labels==0)], 1, mean)
class2_means=apply(expdatafilt_sub_nonlog2[,which(class_labels==1)], 1, mean)
foldchanges=class1_means/class2_means
foldchanges[which(foldchanges<1)]=-1/foldchanges[which(foldchanges<1)]

#Calculate correlations between PRKACA and all other genes
calc_cor_pvalue=function(x,y,method){
	test=cor.test(x=x,y=y,method=method)
	return(test$p.value)
}
expdatafiltmat=as.matrix(expdatafilt)
cor_rvalues_202801_at=apply(expdatafiltmat,1,cor,y=expdatafiltmat["202801_at",],method="spearman")
cor_pvalues_202801_at=apply(expdatafiltmat,1,calc_cor_pvalue,y=expdatafiltmat["202801_at",],method="spearman")
cor_rvalues_216234_s_at=apply(expdatafiltmat,1,cor,y=expdatafiltmat["216234_s_at",],method="spearman")
cor_pvalues_216234_s_at=apply(expdatafiltmat,1,calc_cor_pvalue,y=expdatafiltmat["216234_s_at",],method="spearman")
max_cor_rvalues=apply(cbind(cor_rvalues_202801_at,cor_rvalues_216234_s_at),1,max)
min_cor_pvalues=apply(cbind(cor_pvalues_202801_at, cor_pvalues_216234_s_at),1,min)

#Apply multtest correction to cor.test values
pvalues_202801_at=as.numeric(cor_pvalues_202801_at)
pvalues_adj_202801_at=mt.rawp2adjp(pvalues_202801_at, proc=c("BH"))
pvalues_adj_orig_order_202801_at=pvalues_adj_202801_at$adjp[order(pvalues_adj_202801_at$index),"BH"]
pvalues_216234_s_at=as.numeric(cor_pvalues_216234_s_at)
pvalues_adj_216234_s_at=mt.rawp2adjp(pvalues_216234_s_at, proc=c("BH"))
pvalues_adj_orig_order_216234_s_at=pvalues_adj_216234_s_at$adjp[order(pvalues_adj_216234_s_at$index),"BH"]
min_cor_qvalues=apply(cbind(pvalues_adj_orig_order_202801_at, pvalues_adj_orig_order_216234_s_at),1,min)

#Output results
MTP_summary=cbind(mappingdata[rownames(expdatafilt_sub),],class1_means,class2_means,foldchanges,MTP_results@rawp,MTP_results@adjp,cor_rvalues_202801_at,cor_pvalues_202801_at,pvalues_adj_orig_order_202801_at,cor_rvalues_216234_s_at,cor_pvalues_216234_s_at,pvalues_adj_orig_order_216234_s_at,max_cor_rvalues,min_cor_pvalues,min_cor_qvalues)
colnames(MTP_summary)=c(colnames(mappingdata), "mean (pureFL)", "mean (HCC/Normal)", "fold_change","t_rawp", "t_adjp","rho_vs_202801_at","rho_p_vs_202801_at","rho_q_vs_202801_at","rho_vs_216234_s_at","rho_p_vs_216234_s_at","rho_q_vs_216234_s_at","max_rho","min_rho_p","min_rho_q")
write.table(MTP_summary,file="FL_vs_HCC_Norm_MTPresults.txt",sep="\t",row.names=FALSE)

#Filter down to just probe sets that meet the following criteria
#abs(max_rho)>0.5
#min_rho_q<0.05
#abs(fold_change)>2
#t_adjp<0.05

PRKACA_deregulated=MTP_summary[which(abs(MTP_summary[,"max_rho"])>0.5 & MTP_summary[,"min_rho_q"]<0.05 & abs(MTP_summary[,"fold_change"])>2 & MTP_summary[,"t_adjp"]<0.05),]
PRKACA_upregulated=MTP_summary[which(MTP_summary[,"max_rho"]>0.5 & MTP_summary[,"min_rho_q"]<0.05 & MTP_summary[,"fold_change"]>2 & MTP_summary[,"t_adjp"]<0.05),]

