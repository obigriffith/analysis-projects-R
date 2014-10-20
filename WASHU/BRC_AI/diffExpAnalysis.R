library(multtest)
library(genefilter)
library("gplots")
library("RColorBrewer")
library("biomaRt")


#Load expression data
#setwd("/Users/ogriffit/Dropbox/WashU/Projects/BRC_AI/Expression_data/")
setwd("/Users/obigriffith/Dropbox/WashU/Projects/BRC_AI/Expression_data/")

datafile="SEQPairsGEXMatrixQuery.csv"
rawdata=read.csv(datafile, header = TRUE, na.strings = "NA", row.names=1)

#Load patient and sample info files
patientdatafile="BRC22_patient_info.txt"
sampledatafile="BRC22_sample_info.txt"
sampledata=read.table(sampledatafile, header=TRUE, sep="\t", row.names=1)
patientdata=read.table(patientdatafile, header=TRUE, sep="\t", row.names=1, na.strings="ND")
#Add columns to patient data
patientdata[,"DiffKi67"]=patientdata[,"SurgeryKi67"]-patientdata[,"BaselineKi67"]
patientdata[which(patientdata[,"SurgeryKi67"]<=0.1),"Response"]="Responder"
patientdata[which(patientdata[,"SurgeryKi67"]>0.1),"Response"]="Non-Responder"

#Separate data from annotations
annotations=(rawdata[,1:13])
data=rawdata[,14:length(colnames(rawdata))]

#Get sample and probe names
samples=colnames(data)
substr(samples,1,1)=""
colnames(data)=samples
probes=rownames(data)

#Apply some kind of filtering - These are ratio data so usual Affy filters will not work well
#There is also the problem of what to do with missing values????
na_fun=function(x){
	100*((length(which(is.na(x))))/(length(x)))
}
vars=apply(data,1,var)
percNA=apply(data,1,na_fun)
#Filter out probes with NA values in more than 20% of samples
#Filter out probes without at least some minimum variance 
filt=which(vars > 0.5 & percNA < 20)
filt_data=data[filt,]
filt_annotations=annotations[filt,]
probes_filt=rownames(filt_data)
GeneIDs_filt=filt_annotations[probes_filt,"EntrezGeneID"]
Symbols_filt=as.vector(filt_annotations[probes_filt,"GeneSymbol"])

#Create lists for sample categories
#Note, arbitrarily chose 16178_S over 16178_S_2nd for analysis
Pre=c("15583_BL","15601_BL","15687_BL","15714_BL","15917_BL","15945_BL","15970_BL","16041_BL","16107_BL","16111_BL","16178_BL","16227_BL","16252_BL","16300_BL","16314_BL","16319_BL","16357_BL","16408_BL","16431_BL","16454_BL","16481_BL","16589_BL")
Post=c("15583_S", "15601_S", "15687_S", "15714_S", "15917_S", "15945_S", "15970_S", "16041_S", "16107_S", "16111_S", "16178_S", "16227_S", "16252_S", "16300_S", "16314_S", "16319_S", "16357_S", "16408_S", "16431_S", "16454_S", "16481_S", "16589_S")

#COMPARISON #1 Pre vs Post
samples_sub=c(Pre,Post)
filt_data_sub=filt_data[,samples_sub]
class_labels=c(rep(0,length(Pre)),rep(1,length(Post)))
#Then perform differential expression statistics and multiple testing correction for each comparison
#Note, for paired t-tests, the arrangement of group indices does not matter, as long as the columns are arranged in the same corresponding order between groups. For example, if group 1 is coded as 0, and group 2 is coded as 1, for 3 pairs of data, it does not matter if the label Y is coded as "0,0,0,1,1,1", "1,1,1,0,0,0" "0,1,0,1,0,1" or "1,0,1,0,1,0", the paired differences between groups will be calculated as "group 2 - group 1". See references for more detail.
#Note to get non-parametric equivalent of t-test, use robust=TRUE - not working
MTP_results=MTP(X=filt_data_sub, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.pair", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP", robust=FALSE)
#MTP_results=MTP(X=filt_data_sub, Y=class_labels, na.rm=TRUE, alternative="two.sided", test="t.pair", typeone="fdr", fdr.method="conservative", B=100, method="sd.minP", robust=FALSE)

#To perform standard MU test
var_names_filt=rownames(filt_data_sub)
MU_test_results = array(0, dimnames = list(var_names_filt, c("rawp","adjp")), dim=c(length(var_names_filt),2))
MUfun=function(x){
	pre_values=x[which(class_labels==0)]
	post_values=x[which(class_labels==1)]
	wilcox_result=wilcox.test(x=pre_values, y=post_values, alternative="two.sided", paired=TRUE)
	return(wilcox_result$p.value)	
}
MU_pvalues=apply(filt_data_sub,1,MUfun)

#Perform simple multiple testing correction
pvalues=as.numeric(MU_pvalues)
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
colnames(pvalues_adj_orig_order)=c("MU_rawp","MU_BH_p")

#Calculate basic summary stats (means, fold change, etc)
class1_means=apply (filt_data_sub[,Pre], 1, mean)
class2_means=apply (filt_data_sub[,Post], 1, mean)
foldchanges=class2_means/class1_means
foldchanges[which(foldchanges<1)]=-1/foldchanges[which(foldchanges<1)]
filt_data_sub_log2diff=filt_data_sub[,Post]-filt_data_sub[,Pre]
meanlog2diff=apply(filt_data_sub_log2diff,1,mean)

#Output results.
MTP_summary=cbind(filt_annotations[probes_filt,],class1_means,class2_means,foldchanges,meanlog2diff,MTP_results@rawp, MTP_results@adjp,pvalues_adj_orig_order)
colnames(MTP_summary)=c(colnames(filt_annotations), "mean (Pre)", "mean (Post)", "fold_change", "meanlog2diff","t_rawp", "t_adjp","MU_rawp","MU_BH_p")

#For all 'significant' probes
siggene_summary=MTP_summary[MTP_summary[,"t_adjp"]<0.05,]
siggene_summary2=MTP_summary[MTP_summary[,"MU_BH_p"]<0.05,]

filt_data_sub_sig=filt_data_sub[rownames(siggene_summary),]
filt_data_sub_log2diff_sig=filt_data_sub_log2diff[rownames(siggene_summary),]

#Write final results to files.
write.table(MTP_summary, file="Post_vs_Pre.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(siggene_summary, file="Post_vs_Pre_T_FDR_siggenes.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(siggene_summary2, file="Post_vs_Pre_MU_BH_siggenes.txt", sep="\t", row.names=FALSE, quote=FALSE)


#COMPARISON #2 Responders vs non-responders
#Apply filtering to just post-surgery subset of data
samples_sub2=c(Post)
data_post=data[,samples_sub2]

vars2=apply(data_post,1,var)
percNA2=apply(data_post,1,na_fun)
#Filter out probes with NA values in more than 20% of samples
#Filter out probes without at least some minimum variance 
filt2=which(vars2 > 0.5 & percNA2 < 20)
filt_data2=data_post[filt2,]
filt_annotations2=annotations[filt2,]
probes_filt2=rownames(filt_data2)
GeneIDs_filt2=filt_annotations2[probes_filt2,"EntrezGeneID"]
Symbols_filt2=as.vector(filt_annotations2[probes_filt2,"GeneSymbol"])
filt_data_sub2=filt_data2[,samples_sub2]

class_labels2=patientdata[sampledata[samples_sub2,"CommonName"],"Response"]
class_labels2[which(class_labels2=="Non-Responder")]=0
class_labels2[which(class_labels2=="Responder")]=1

#To perform standard MU test
var_names_filt2=rownames(filt_data_sub2)
MU_test_results = array(0, dimnames = list(var_names_filt2, c("rawp","adjp")), dim=c(length(var_names_filt2),2))
MUfun=function(x){
	nonresponder_values=x[which(class_labels2==0)]
	responder_values=x[which(class_labels2==1)]
	wilcox_result=wilcox.test(x=nonresponder_values, y=responder_values, alternative="two.sided", paired=FALSE)
	return(wilcox_result$p.value)	
}
MU_pvalues2=apply(filt_data_sub2,1,MUfun)

#Perform simple multiple testing correction
pvalues2=as.numeric(MU_pvalues2)
pvalues_adj2=mt.rawp2adjp(pvalues2, proc=c("BH"))
pvalues_adj_orig_order2=pvalues_adj2$adjp[order(pvalues_adj2$index),]
colnames(pvalues_adj_orig_order2)=c("MU_rawp","MU_BH_p")

MTP_results2=MTP(X=filt_data_sub2, Y=class_labels2, na.rm=TRUE, alternative="two.sided", test="t.twosamp.equalvar", typeone="fdr", fdr.method="conservative", B=1000, method="sd.minP", robust=TRUE)

#Calculate basic summary stats (means, fold change, etc)
class1_means2=apply (filt_data_sub2[,which(class_labels2==0)], 1, mean)
class2_means2=apply (filt_data_sub2[,which(class_labels2==1)], 1, mean)
foldchanges2=class2_means2/class1_means2
foldchanges2[which(foldchanges2<1)]=-1/foldchanges2[which(foldchanges2<1)]
meanlog2diff2=class2_means2-class1_means2

#Output results.
MTP_summary2=cbind(filt_annotations2[probes_filt2,],class1_means2,class2_means2,foldchanges2,meanlog2diff2,MTP_results2@rawp, MTP_results2@adjp,pvalues_adj_orig_order2)
colnames(MTP_summary2)=c(colnames(filt_annotations2), "mean (Non-Responder)", "mean (Responder)", "fold_change", "meanlog2diff","MU_MTP_rawp", "MU_MTP_adjp","MU_rawp","MU_BH_p")

#For all 'significant' probes
siggene_summary3=MTP_summary2[MTP_summary2[,"MU_rawp"]<0.05,]
siggene_summary4=MTP_summary2[MTP_summary2[,"MU_MTP_adjp"]<0.05,]
filt_data_sub_sig2=filt_data_sub2[rownames(siggene_summary3),]

#Write final results to files.
write.table(MTP_summary2, file="Responder_vs_NonResponder.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(siggene_summary3, file="Responder_vs_NonResponder_siggenes.txt", sep="\t", row.names=FALSE, quote=FALSE)


###Heatmaps###
source("/Users/obigriffith/Dropbox/drug_predictors/Rscripts/heatmap.3.R")

#Heatmap #1: All genes passing basic filter
sample_common_names=sampledata[samples_sub,"CommonName"]

#Create color side bars
colors=rainbow(6)
subtype_colors=as.vector(patientdata[sample_common_names,"BaselineSubtype"])
subtype_colors[which(subtype_colors=="LumA")]=colors[1]
subtype_colors[which(subtype_colors=="LumB")]=colors[2]

PrePost_colors=as.vector(sampledata[samples_sub,"Pre_Post"])
PrePost_colors[which(PrePost_colors=="Pre")]=colors[3]
PrePost_colors[which(PrePost_colors=="Post")]=colors[4]

Responder_colors=as.vector(patientdata[sample_common_names,"Response"])
Responder_colors[which(Responder_colors=="Responder")]=colors[5]
Responder_colors[which(Responder_colors=="Non-Responder")]=colors[6]

Ki67_colors=c(patientdata[sampledata[Pre,"CommonName"],"BaselineKi67"],patientdata[sampledata[Post,"CommonName"],"SurgeryKi67"])
Ki67_colors=(round(Ki67_colors, digits=2)*100)+1
blues=colorRampPalette(brewer.pal(9,"Blues"))(101)
Ki67_colors=blues[Ki67_colors]
Ki67_colors[which(is.na(Ki67_colors))]="grey"

patient_colors=rainbow(22)
names(patient_colors)=levels(sample_common_names)
sample_patient_colors=patient_colors[sample_common_names]

clab=cbind(sample_patient_colors,Ki67_colors,Responder_colors,PrePost_colors,subtype_colors)
colnames(clab)=c("Patient","Ki67","Response","Pre/Post","Subtype")

#Create a custom color palatte for heatmap from yellow (down) through white (no diff) to blue (up)
nHalf=50
Min=min(filt_data_sub)
Max=max(filt_data_sub)
Thresh=0
rc1 <- colorRampPalette(colors = c("yellow", "white"), space="Lab")(nHalf)    
rc2 <- colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

pdf(file="BRC22_all_filtered_probes_heatmap.pdf")
main_title=paste("BRC22 pre/post, all",length(rownames(filt_data_sub)),"probes passing filters")
par(cex.main=0.9)
heatmap.3(as.matrix(filt_data_sub), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(4,10), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=FALSE, labCol=sample_common_names, cexRow=1, col=rampcols, breaks=rampbreaks, NumColSideColors=5, KeyValueName="Log2 ratio")
legend("topright",legend=c("LumA","LumB","","Pre-Treat","Post-Treat","","Responder","Non-Responder","","Ki67 (%)","High (100)","   .","   .","   .","Low (0)","N/A"), fill=c(colors[1],colors[2],"white",colors[3],colors[4],"white",colors[5],colors[6],"white","white",blues[100],blues[75],blues[50],blues[25],blues[1],"grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.5)
dev.off()


#Heatmap #2: Differences between paired Pre/Post samples for sig genes
samples_sub_diff=colnames(filt_data_sub_log2diff_sig)
patient_common_names=sampledata[samples_sub_diff,"CommonName"]

#Create color side bars
colors=rainbow(6)
subtype_colors=as.vector(patientdata[patient_common_names,"BaselineSubtype"])
subtype_colors[which(subtype_colors=="LumA")]=colors[1]
subtype_colors[which(subtype_colors=="LumB")]=colors[2]

Responder_colors=as.vector(patientdata[patient_common_names,"Response"])
Responder_colors[which(Responder_colors=="Responder")]=colors[5]
Responder_colors[which(Responder_colors=="Non-Responder")]=colors[6]

Ki67_colors=patientdata[sampledata[samples_sub_diff,"CommonName"],"SurgeryKi67"]
Ki67_colors=(round(Ki67_colors, digits=2)*100)+1
blues=colorRampPalette(brewer.pal(9,"Blues"))(101)
Ki67_colors=blues[Ki67_colors]
Ki67_colors[which(is.na(Ki67_colors))]="grey"

clab=cbind(Ki67_colors,Responder_colors,subtype_colors)
colnames(clab)=c("Ki67","Response","Subtype")

#Create a custom color palatte for heatmap from yellow (down) through white (no diff) to blue (up)
nHalf=50
Min=min(filt_data_sub_log2diff_sig)
Max=max(filt_data_sub_log2diff_sig)
Thresh=0
rc1 <- colorRampPalette(colors = c("yellow", "white"), space="Lab")(nHalf)    
rc2 <- colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#Filter down to just probes with a GeneSymbol
filt_data_sub_log2diff_sig_wSymb=filt_data_sub_log2diff_sig[which(annotations[rownames(filt_data_sub_log2diff_sig),"GeneSymbol"]!=""),]
gene_names=as.vector(annotations[rownames(filt_data_sub_log2diff_sig_wSymb),"GeneSymbol"])

pdf(file="BRC22_all_filtered_sig_probes_diff_heatmap.pdf")
main_title=paste("BRC22 post-pre diff,",length(rownames(filt_data_sub_log2diff_sig_wSymb)),"most significant probes (genes)")
par(cex.main=0.8)
heatmap.3(as.matrix(filt_data_sub_log2diff_sig_wSymb), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(4,10), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=patient_common_names, labRow=gene_names, cexRow=0.4, col=rampcols, breaks=rampbreaks, NumColSideColors=3, KeyValueName="Diff(Log2 ratio)",)
legend("topright",legend=c("LumA","LumB","","Responder","Non-Responder","","Ki67 (%)","High (100)","   .","   .","   .","Low (0)"), fill=c(colors[1],colors[2],"white",colors[5],colors[6],"white","white",blues[100],blues[75],blues[50],blues[25],blues[1]), border=FALSE, bty="n", y.intersp = 0.7, cex=0.5)
dev.off()

#Heatmap #3: Differences between responders and non-responders for sig genes
patient_common_names=sampledata[samples_sub2,"CommonName"]

#Create color side bars
colors=rainbow(6)
subtype_colors=as.vector(patientdata[patient_common_names,"BaselineSubtype"])
subtype_colors[which(subtype_colors=="LumA")]=colors[1]
subtype_colors[which(subtype_colors=="LumB")]=colors[2]

Responder_colors=as.vector(patientdata[patient_common_names,"Response"])
Responder_colors[which(Responder_colors=="Responder")]=colors[5]
Responder_colors[which(Responder_colors=="Non-Responder")]=colors[6]

Ki67_colors=patientdata[sampledata[samples_sub2,"CommonName"],"SurgeryKi67"]
Ki67_colors=(round(Ki67_colors, digits=2)*100)+1
blues=colorRampPalette(brewer.pal(9,"Blues"))(101)
Ki67_colors=blues[Ki67_colors]
Ki67_colors[which(is.na(Ki67_colors))]="grey"

clab=cbind(Ki67_colors,Responder_colors,subtype_colors)
colnames(clab)=c("Ki67","Response","Subtype")

#Create a custom color palatte for heatmap from yellow (down) through white (no diff) to blue (up)
nHalf=50
Min=min(filt_data_sub_sig2)
Max=max(filt_data_sub_sig2)
Thresh=0
rc1 <- colorRampPalette(colors = c("yellow", "white"), space="Lab")(nHalf)    
rc2 <- colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#Filter down to just probes with a GeneSymbol
filt_data_sub_sig2_wSymb= filt_data_sub_sig2[which(annotations[rownames(filt_data_sub_sig2),"GeneSymbol"]!=""),]
gene_names=as.vector(annotations[rownames(filt_data_sub_sig2_wSymb),"GeneSymbol"])

pdf(file="BRC22_all_filtered_sig_probes_Resp_vs_NonResp_heatmap.pdf")
main_title=paste("BRC22 Resp vs NonResp,",length(rownames(filt_data_sub_sig2_wSymb)),"most sig. probes (genes)")
par(cex.main=0.8)
heatmap.3(as.matrix(filt_data_sub_sig2_wSymb), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(4,10), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labCol=patient_common_names, labRow=gene_names, cexRow=0.4, col=rampcols, breaks=rampbreaks, NumColSideColors=3, KeyValueName="Log2 ratio",)
legend("topright",legend=c("LumA","LumB","","Responder","Non-Responder","","Ki67 (%)","High (100)","   .","   .","   .","Low (0)"), fill=c(colors[1],colors[2],"white",colors[5],colors[6],"white","white",blues[100],blues[75],blues[50],blues[25],blues[1]), border=FALSE, bty="n", y.intersp = 0.7, cex=0.5)
dev.off()








