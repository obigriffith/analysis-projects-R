library("gplots")

#Load variant data
coverage_matrix_file="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/bam_readcounts_all_bams/final_coverage_matrix.tsv"
refcount_matrix_file="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/bam_readcounts_all_bams/final_refcount_matrix.tsv"
vaf_matrix_file="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/bam_readcounts_all_bams/final_vaf_matrix.tsv"
varcount_matrix_file="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/bam_readcounts_all_bams/final_varcount_matrix.tsv"
coverage_matrix_data=read.table(coverage_matrix_file,sep="\t",header=TRUE, as.is=c(1:7), check.names=FALSE)
refcount_matrix_data=read.table(refcount_matrix_file,sep="\t",header=TRUE, as.is=c(1:7), check.names=FALSE)
vaf_matrix_data=read.table(vaf_matrix_file,sep="\t",header=TRUE, as.is=c(1:7), check.names=FALSE)
varcount_matrix_data=read.table(varcount_matrix_file,sep="\t",header=TRUE, as.is=c(1:7), check.names=FALSE)

#Load sample info
sample_info_file="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/sample_info.txt"
sample_info_data=read.table(sample_info_file,sep="\t",header=TRUE, as.is=c(1:4))
sample_ids=sample_info_data[,"SampleID"]
sample_names=sample_info_data[,"SampleName"]

#Set result dir
result_dir="/Users/ogriffit/Dropbox/WashU/Projects/BRAF/bam_readcounts_all_bams/"
setwd(result_dir)

#Set outfiles
VAFheatmapfile="BRAF_all_lines_RNA_vs_DNA_VAF.pdf"
VAFDNAheatmapfile="BRAF_all_lines_DNA_VAF.pdf"
VAFRNAheatmapfile="BRAF_all_lines_RNA_VAF.pdf"
VAFfilteredheatmapfile="BRAF_all_lines_RNA_vs_DNA_VAF_filt.pdf"
VAFstringentfilteredheatmapfile="BRAF_all_lines_RNA_vs_DNA_VAF_filt_stringent.pdf"
VAFDNAstringentfilteredheatmapfile="BRAF_all_lines_DNA_VAF_filt_stringent.pdf"

#Create Heatmaps
#Use modified heatmap.2 command to allow multiple color side bars
source("/Users/ogriffit/Dropbox/drug_predictors/Rscripts/heatmap.3.R")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap #1 - DNA + RNA, VAF
#Extract just data columns for clustering purposes
vafdata=as.matrix(vaf_matrix_data[,sample_ids])
coveragedata=as.matrix(coverage_matrix_data[,sample_ids])

#Set up colors for sidebars
RNADNAcolors=sample_info_data[,"RNA_DNA"]
RNADNAcolors[RNADNAcolors=="RNA"]="red"
RNADNAcolors[RNADNAcolors=="DNA"]="blue"
Samplecolors=c(rainbow(10),rainbow(10))
clab=cbind(RNADNAcolors,Samplecolors)
colnames(clab)=c("RNA/DNA","Sample")

pdf(file=VAFheatmapfile, height=10, width=7.5)
main_title="BRAF cell lines, VAF, RNA vs DNA"
par(cex.main=1)
heatmap.3(vafdata, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(7,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
legend("left",legend=c("RNA","DNA"), fill=c("red","blue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#Heatmap #2 - DNA only, VAF
#Extract just data columns for clustering purposes
DNA_sample_ids=sample_info_data[which(sample_info_data[,"RNA_DNA"]=="DNA"),"SampleID"]
DNAvafdata=as.matrix(vaf_matrix_data[,DNA_sample_ids])

#Set up colors for sidebars
Samplecolors=rainbow(10)
clab=cbind(Samplecolors,Samplecolors)
colnames(clab)=c("Sample","")

pdf(file=VAFDNAheatmapfile, height=10, width=7.5)
main_title="BRAF cell lines, VAF, DNA"
par(cex.main=1)
heatmap.3(DNAvafdata, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(9,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
dev.off()

#Heatmap #3 - RNA only, VAF, only variants with minimum coverage/VAF in at least one RNA sample
#Extract just data columns for clustering purposes
RNA_sample_ids=sample_info_data[which(sample_info_data[,"RNA_DNA"]=="RNA"),"SampleID"]
RNAvafdata=as.matrix(vaf_matrix_data[,RNA_sample_ids])
RNAcoveragedata=as.matrix(coverage_matrix_data[,RNA_sample_ids])

#Define a function to count number of values passing true test
count_true_fun=function(x){
	count=length(which(x))
	return(count)
}

RNA_filt_test=RNAvafdata>10 & RNAcoveragedata>10
RNA_filt_count_data=apply(RNA_filt_test,1,count_true_fun)
RNA_filt_passed=which(RNA_filt_count_data>=1)
RNAvafdata_filt=RNAvafdata[RNA_filt_passed,]

#Set up colors for sidebars
Samplecolors=rainbow(10)
clab=cbind(Samplecolors,Samplecolors)
colnames(clab)=c("Sample","")

pdf(file=VAFRNAheatmapfile, height=10, width=7.5)
main_title="BRAF cell lines, VAF, RNA, min 10% VAF & 10X"
par(cex.main=1)
heatmap.3(RNAvafdata_filt, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(9,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
dev.off()

#Heatmap #4 - DNA + RNA, VAF, only variants with minimum coverage/VAF in at least one RNA sample
#Extract just data columns for clustering purposes
vafdata_filt=as.matrix(vaf_matrix_data[RNA_filt_passed,sample_ids])

#Set up colors for sidebars
Samplecolors=c(rainbow(10),rainbow(10))
clab=cbind(RNADNAcolors,Samplecolors)
colnames(clab)=c("RNA/DNA","Sample")

pdf(file=VAFfilteredheatmapfile, height=10, width=7.5)
main_title="BRAF cell lines, VAF, RNA vs DNA, min RNA 10% VAF & 10X"
par(cex.main=1)
heatmap.3(vafdata_filt, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(9,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
legend("left",legend=c("RNA","DNA"), fill=c("red","blue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#Heatmap #5
#More stringent filtering strategy for DNA vs RNA analysis
#1 - eliminate any VAF estimates based on coverage less than 10X for RNA or DNA (set to 0)
vafdata_filt=as.matrix(vaf_matrix_data[,sample_ids])
coveragedata=as.matrix(coverage_matrix_data[,sample_ids])
vafdata_filt[which(coveragedata<10,arr.ind=TRUE)]=0

#2 - eliminate any variants with VAF > 5% in parental lines
parental_VAF_test=vafdata_filt[,c(1,7,11,17)]>5
parental_VAF_test_count_data=apply(parental_VAF_test,1,count_true_fun)
parental_VAF_test_passed=which(parental_VAF_test_count_data==0)
vafdata_filt=vafdata_filt[parental_VAF_test_passed,]

#3 - eliminate any VAF estimates below 10% or above 90%
vafdata_filt[which((vafdata_filt<10 | vafdata_filt>90), arr.ind=TRUE)]=0

#4 - eliminate any variants without at least one RNA sample having VAF >= 10%
RNA_filt_test=vafdata_filt[,RNA_sample_ids]>=10
RNA_filt_count_data=apply(RNA_filt_test,1,count_true_fun)
RNA_filt_passed2=which(RNA_filt_count_data>=1)
vafdata_filt=vafdata_filt[RNA_filt_passed2,]

pdf(file=VAFstringentfilteredheatmapfile, height=10, width=7.5)
main_title="VAF, RNA vs DNA, 10x, non-parent, 10<VAF<90, RNA VAF>10"
par(cex.main=1)
heatmap.3(vafdata_filt, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(9,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
legend("left",legend=c("RNA","DNA"), fill=c("red","blue"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()



#Heatmap #6
#More stringent filtering strategy for DNA only analysis
#1 - eliminate any VAF estimates based on coverage less than 10X for DNA (set to 0)
vafdata_filt=as.matrix(vaf_matrix_data[,DNA_sample_ids])
coveragedata=as.matrix(coverage_matrix_data[,DNA_sample_ids])
vafdata_filt[which(coveragedata<10,arr.ind=TRUE)]=0

#2 - eliminate any variants with VAF > 5% in parental lines
parental_VAF_test=vafdata_filt[,c(1,7)]>5
parental_VAF_test_count_data=apply(parental_VAF_test,1,count_true_fun)
parental_VAF_test_passed=which(parental_VAF_test_count_data==0)
vafdata_filt=vafdata_filt[parental_VAF_test_passed,]

#3 - eliminate any VAF estimates below 10% or above 90%
vafdata_filt[which((vafdata_filt<10 | vafdata_filt>90), arr.ind=TRUE)]=0

#4 - eliminate any variants without at least one DNA sample having VAF >= 10%
DNA_filt_test=vafdata_filt[,DNA_sample_ids]>=10
DNA_filt_count_data=apply(DNA_filt_test,1,count_true_fun)
DNA_filt_passed2=which(DNA_filt_count_data>=1)
vafdata_filt=vafdata_filt[DNA_filt_passed2,]

#Set up colors for sidebars
Samplecolors=sample_info_data[1:10,"Group"]
Samplecolors[Samplecolors=="VacoParent"]="lightblue"
Samplecolors[Samplecolors=="Vaco"]="blue"
Samplecolors[Samplecolors=="F6Parent"]="lightgreen"
Samplecolors[Samplecolors=="F6"]="green"
clab=cbind(Samplecolors,Samplecolors)
colnames(clab)=c("","")

pdf(file=VAFDNAstringentfilteredheatmapfile, height=10, width=10)
main_title="VAF, DNA, 10x, non-parent, 10<VAF<90, DNA 1 VAF>10"
par(cex.main=1)
heatmap.3(vafdata_filt, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="column", margins=c(9,4), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=sample_names, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)), NumColSideColors=2, KeyValueName="VAF")
dev.off()







