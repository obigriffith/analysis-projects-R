library("heatmap.plus")
library(genefilter)
library(gplots)
library(mclust)

#datafile="C:/Users/Obi/Documents/Projects/Cepheid/processing/processed/customCDF/ALL_gcrma.txt"
#datafile="C:/Users/Obi/Documents/Projects/Cepheid/processing/processed/standardCDF/ALL_gcrma.txt"
#datafile="C:/Users/Obi/Documents/Projects/Cepheid/processing/processed2/customCDF/ALL_gcrma.txt"
datafile="C:/Users/Obi/Documents/Projects/Cepheid/processing/processed2/standardCDF/ALL_gcrma.txt"

#clindatafile="C:/Users/Obi/Documents/Projects/Cepheid/clinical_data/filtered_combined_data_anno.txt"
clindatafile="C:/Users/Obi/Documents/Projects/Cepheid/clinical_data/filtered_combined_data_anno.3.txt"

#outdir="C:/Users/Obi/Documents/Projects/Cepheid/analyzing/analysis/customCDF/ALL_samples"
#outdir="C:/Users/Obi/Documents/Projects/Cepheid/analyzing/analysis/standardCDF/ALL_samples"
#outdir="C:/Users/Obi/Documents/Projects/Cepheid/analyzing/analysis2/customCDF/ALL_samples"
outdir="C:/Users/Obi/Documents/Projects/Cepheid/analyzing/analysis2/standardCDF/ALL_samples"

#Results files
duplicates_file = "ALL_gcrma_duplicate_samples.txt"
duplicates_pdf = "ALL_gcrma_sample_correlations.pdf"
outfile = "ALL_gcrma_heatmap_filtered_genes.pdf"
outfile2 = "ALL_gcrma_heatmap_ERBB2_GRB7_clustered.pdf"
outfile2B = "ALL_gcrma_heatmap_ERBB2_GRB7_ESR1_clustered.pdf"
outfile2C = "ALL_gcrma_heatmap_ERBB2_Amplicon_clustered.pdf"
outfile2D = "ALL_gcrma_heatmap_ERBB2_Amplicon_clean_clustered.pdf"
outfile2E = "ALL_gcrma_heatmap_ERBB2_Amplicon_clustered_chr_order.pdf"
outfile3 = "ALL_gcrma_heatmap_ERBB2_ordered.pdf"
outfile3B = "ALL_gcrma_heatmap_ESR1_ordered.pdf"
outfile3C = "ALL_gcrma_heatmap_ESR1_allprobes.pdf"
outfile3D = "ALL_gcrma_heatmap_ERBB2_allprobes.pdf"
outfile4 = "ALL_gcrma_waterfall_ERBB2.pdf"
outfile4B = "ALL_gcrma_waterfall_ESR1.pdf"
outfile5A = "ALL_gcrma_hist_ERBB2_mix_model.pdf"
outfile5B = "ALL_gcrma_hist_ERBB2.pdf"
outfile5C = "ALL_gcrma_hist_ESR1_mix_model.pdf"
outfile5D = "ALL_gcrma_hist_ESR1.pdf"
outfile5E = "ALL_gcrma_hist_ERBB2amp_mix_model.pdf"
outfile5F = "ALL_gcrma_hist_ERBB2amp.pdf"
outfile6A = "ALL_gcrma_ERBB2_pos_samples.txt"
outfile6B = "ALL_gcrma_ERBB2amp_pos_samples.txt"
outfile6C = "ALL_gcrma_ESR1_neg_samples.txt"
newclinfile = "filtered_combined_data_anno.3.array_annots.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clin_data=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
raw_data=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

#Change to output dir
setwd(outdir)
header=colnames(raw_data)

#Preliminary gene filtering
X=raw_data[,4:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=raw_data[filt,] 

#Use this filtered subset of the data to look for duplicates in the data
x=as.matrix(filt_Data[,4:length(header)])
#x=as.matrix(raw_data[,4:length(header)]) #Checked to make sure using entire probeset produces same results. It did.
cors=cor(x, method="pearson")
cors[lower.tri(cors, diag=TRUE)]=NA #Set one triangle, plus the diagonal to NA
pdf(file=duplicates_pdf)
hist(cors, col="blue", breaks=1000, xlim=c(-0.2,1), main="Sample correlations for all 998 x 998 samples, variant probes", xlab="Pearson cor")
dev.off()
dup_count1=length(which(cors>0.99)) #Check how many pairs have correlation of 1 (use >0.99, otherwise miss cases that are 0.99999999)
dup_count1

duplicate_pairs=matrix(data = NA, nrow = dup_count1, ncol = 2)

#Go through and identify sample pairs with cor==1
dup_count2=0
samples=rownames(cors)
for (i in 1:length(samples)){
  sample=samples[i]
  dup_test=which(cors[i,]>0.99)
  if (length(dup_test)>0){
    dup_count2=dup_count2+1
    dup_samples=names(dup_test)
    print(paste(i,sample,dup_samples))
    duplicate_pairs[dup_count2,1]=sample
    duplicate_pairs[dup_count2,2]=dup_samples
  }
}

dup_count2 #Make sure you found all the duplicates
write.table(duplicate_pairs, file=duplicates_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#Single-gene Data
#ER data
ESR1_Data_all=raw_data[which(raw_data[,3]=="ESR1"),]

#ERBB2 Amplicon genes from Kauraniemi and Kallioniemi 2006.
ERBB2_Data_all=raw_data[which(raw_data[,3]=="ERBB2"),]
GRB7_Data=raw_data[which(raw_data[,3]=="GRB7"),]
NEUROD2_Data=raw_data[which(raw_data[,3]=="NEUROD2"),]
#PPP1R1B_Data=raw_data[which(raw_data[,3]=="PPP1R1B"),] #No probes found
STARD3_Data=raw_data[which(raw_data[,3]=="STARD3"),]
TCAP_Data=raw_data[which(raw_data[,3]=="TCAP"),]
PNMT_Data=raw_data[which(raw_data[,3]=="PNMT"),]
PGAP3_Data=raw_data[which(raw_data[,3]=="PGAP3"),] #AKA PERLD1
#C17orf37_Data=raw_data[which(raw_data[,3]=="C17orf37"),] #No probes found
IKZF3_Data=raw_data[which(raw_data[,3]=="IKZF3"),] #AKA ZNFN1A3

#Grouped 
ERBB2_GRB7_Data=rbind(ERBB2_Data_all,GRB7_Data)
ERBB2_GRB7_ESR1_Data=rbind(ERBB2_Data_all,GRB7_Data,ESR1_Data_all)
ERBB2_Amplicon_Data=rbind(NEUROD2_Data,STARD3_Data,TCAP_Data,PNMT_Data,PGAP3_Data,ERBB2_Data_all,GRB7_Data,IKZF3_Data)


#Replace gene-level data with best probe (if multiple probes) (for standardCDF)
ERBB2_Data=raw_data[which(raw_data[,1]=="216836_s_at"),] 
ESR1_Data=raw_data[which(raw_data[,1]=="205225_at"),] 
PGAP3_Data=raw_data[which(raw_data[,1]=="55616_at"),]

#Filter down to just reliable probes for ERBB2 amplicon
ERBB2_Amplicon_Data_clean=rbind(ERBB2_Data,GRB7_Data,STARD3_Data,PGAP3_Data)
ERBB2_Amplicon_Data_clean2=rbind(NEUROD2_Data,STARD3_Data,TCAP_Data,PNMT_Data,PGAP3_Data,ERBB2_Data,GRB7_Data,IKZF3_Data)

#Set probe_name for ERBB2 and ESR1 to use in plot titles
ERBB2_probe=ERBB2_Data[1,1]
ESR1_probe=ESR1_Data[1,1]

#Set up colors for sidebars
colors=rainbow(8)
#Study
study_colors=as.vector(clin_data[,1])
study_colors[study_colors=="Desmedt_2007"]=colors[1]
study_colors[study_colors=="Upsalla_combined"]=colors[2]
study_colors[study_colors=="Loi_2007"]=colors[3]
study_colors[study_colors=="Schmidt_2008"]=colors[4]
study_colors[study_colors=="Sotiriou_2006"]=colors[5]
study_colors[study_colors=="Symmans_2010"]=colors[6]
study_colors[study_colors=="Wang_2005"]=colors[7]
study_colors[study_colors=="Zhang_2009"]=colors[8]

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}


#Heatmap (single color sidebar) - All filtered genes
x=as.matrix(filt_Data[,4:length(header)])
pdf(file=outfile)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=study_colors, labRow=FALSE, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ERBB2/GRB7 only
x=as.matrix(ERBB2_GRB7_Data[,4:length(header)])
probe_names=ERBB2_GRB7_Data[,1]
gene_names=ERBB2_GRB7_Data[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile2)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="column", Rowv=FALSE, ColSideColors=study_colors, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ERBB2/GRB7/ESR1 only
x=as.matrix(ERBB2_GRB7_ESR1_Data[,4:length(header)])
probe_names=ERBB2_GRB7_ESR1_Data[,1]
gene_names=ERBB2_GRB7_ESR1_Data[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile2B)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="column", Rowv=FALSE, ColSideColors=study_colors, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ERBB2 Amplicon genes only
x=as.matrix(ERBB2_Amplicon_Data[,4:length(header)])
probe_names=ERBB2_Amplicon_Data[,1]
gene_names=ERBB2_Amplicon_Data[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile2C)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=study_colors, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()




#Heatmap (single color sidebar) - ERBB2 Amplicon genes only, ordered by chromosome, best probe for each gene only
x=as.matrix(ERBB2_Amplicon_Data_clean2[,4:length(header)])
probe_names=ERBB2_Amplicon_Data_clean2[,1]
gene_names=ERBB2_Amplicon_Data_clean2[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile2E)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, dendrogram="column", Rowv=FALSE, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=study_colors, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()







#Heatmap (single color sidebar) - ERBB2 Amplicon genes only, good probes only
x=as.matrix(ERBB2_Amplicon_Data_clean[,4:length(header)])
probe_names=ERBB2_Amplicon_Data_clean[,1]
gene_names=ERBB2_Amplicon_Data_clean[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile2D)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=study_colors, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ERBB2 only, order samples by decreasing ERBB2 instead of clustering. Color side bar requires same reorder
x=sort(as.matrix(ERBB2_Data[,4:length(header)]), decreasing=TRUE)
z=order(as.matrix(ERBB2_Data[,4:length(header)]), decreasing=TRUE)
study_colors_reorder=study_colors[z]
x=rbind(x,x) #duplicate to trick heatmap.2 which requires minimum of 2 x 2 matrix
pdf(file=outfile3)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE, ColSideColors=study_colors_reorder, main=paste("ERBB2 -",ERBB2_probe), labRow=FALSE, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ESR1 only, order samples by decreasing ESR1 instead of clustering. Color side bar requires same reorder
x=sort(as.matrix(ESR1_Data[,4:length(header)]), decreasing=TRUE)
z=order(as.matrix(ESR1_Data[,4:length(header)]), decreasing=TRUE)
study_colors_reorder=study_colors[z]
x=rbind(x,x) #duplicate to trick heatmap.2 which requires minimum of 2 x 2 matrix
pdf(file=outfile3B)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE, ColSideColors=study_colors_reorder, main=paste("ESR1 -",ESR1_probe), labRow=FALSE, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Heatmap (single color sidebar) - ESR1 only, all probes
x=as.matrix(ESR1_Data_all[,4:length(header)])
probe_names=ESR1_Data_all[,1]
gene_names=ESR1_Data_all[,3]
row_names=paste(probe_names,"(",gene_names,")",sep="")
pdf(file=outfile3C)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=study_colors, main="ESR1", labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()





#Waterfall plot - ERBB2
y=sort(as.numeric(ERBB2_Data[,4:length(header)]), decreasing=TRUE)
x=0:(length(y)-1)
pdf(file=outfile4)
plot(y, type="n", main=paste("ERBB2 (",ERBB2_probe,") waterfall for all 998 samples", sep=""), ylab="Log2 GCRMA value", xlab="Sample")
polygon(c(min(x), x, max(x), 0), c(0, y, 0, 0), col="blue")
dev.off()

#Waterfall plot - ESR1
y=sort(as.numeric(ESR1_Data[,4:length(header)]), decreasing=TRUE)
x=0:(length(y)-1)
pdf(file=outfile4B)
plot(y, type="n", main=paste("ESR1 (",ESR1_probe,") waterfall for all 998 samples", sep=""), ylab="Log2 GCRMA value", xlab="Sample")
polygon(c(min(x), x, max(x), 0), c(0, y, 0, 0), col="blue")
dev.off()

#Attempt to separate ERBB2 data based on gaussian mixture model
#ERBB2 alone
#Allow mclust to determine best model, then filter out extreme dist'n (Want to filter ERBB2+, therefore -> right/last component)
x=as.numeric(ERBB2_Data[,4:length(header)])
mclust_ERBB2=Mclust(x)
summary(mclust_ERBB2, x) #Gives you list of values returned
classification_ERBB2=mclust_ERBB2$classification
num_clusters=mclust_ERBB2$G
pdf(file=outfile5A)
par(mfrow=c(2,1), oma=c(2,2,2,2))
hist(x[classification_ERBB2<num_clusters], xlim=c(2,15), col="blue", xlab="Log2 GCRMA value", main="ERBB2 components to retain")
hist(x[classification_ERBB2==num_clusters], xlim=c(2,15), col="red", xlab="Log2 GCRMA value", main="ERBB2 components to filter")
title(main="Separation of ERBB2+ by model-based clustering", outer=TRUE)
ERBB2_cutoff=min(x[classification_ERBB2==num_clusters]) #Use minimum value of last cluster as cutoff point
ERBB2_filtered_count=length(which(x>=ERBB2_cutoff))
dev.off()

#Histogram - ERBB2
pdf(file=outfile5B)
hist(as.numeric(ERBB2_Data[,4:length(header)]), xlim=c(2,15), breaks=40, col="blue", main=paste("ERBB2 (",ERBB2_probe,") histogram for all 998 samples", sep=""), xlab="Log2 GCRMA value", ylab="Frequency")
abline(v=ERBB2_cutoff, col="red")
legend("left", legend=c(paste("cutoff=",ERBB2_cutoff,sep=""),paste("N=",ERBB2_filtered_count,sep="")),bty="n")
dev.off()

#Print out filtered samples to file
x=as.matrix(ERBB2_Data[,4:length(header)])
samples=colnames(x)
ERBB2_filter_samples=samples[which(x>=ERBB2_cutoff)]
write.table(ERBB2_filter_samples, file=outfile6A, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#ERBB2 Amplicon RankSum
#Allow mclust to determine best model, then filter out extreme dist'n (Want to filter ERBB2+, therefore -> right/last component)
x=as.matrix(ERBB2_Amplicon_Data_clean[,4:length(header)])
x=apply(x,1,rank)
x=apply(x,1,sum)
mclust_ERBB2_RS=Mclust(x)
summary(mclust_ERBB2_RS, x) #Gives you list of values returned
classification_ERBB2_RS=mclust_ERBB2_RS$classification
num_clusters=mclust_ERBB2_RS$G
pdf(file=outfile5E)
par(mfrow=c(2,1), oma=c(2,2,2,2))
hist(x[classification_ERBB2_RS<num_clusters], xlim=c(100,4000), col="blue", xlab="RankSum (4 informative probes)", main="ERBB2_amp components to retain")
hist(x[classification_ERBB2_RS==num_clusters], xlim=c(100,4000), col="red", xlab="RankSum (4 informative probes)", main="ERBB2_amp components to filter")
title(main="Separation of ERBB2_amp+ by model-based clustering", outer=TRUE)
ERBB2_RS_cutoff=min(x[classification_ERBB2_RS==num_clusters]) #Use minimum value of last cluster as cutoff point
ERBB2_RS_filtered_count=length(which(x >= ERBB2_RS_cutoff))
dev.off()

#Histogram - ERBB2 Amplicon RankSum
pdf(file=outfile5F)
hist(as.numeric(x), xlim=c(100,4000), breaks=100, col="blue", main="ERBB2_amp histogram for all 998 samples", xlab="RankSum (4 informative probes)", ylab="Frequency")
abline(v=ERBB2_RS_cutoff, col="red")
legend("left", legend=c(paste("cutoff=",ERBB2_RS_cutoff,sep=""),paste("N=",ERBB2_RS_filtered_count,sep="")),bty="n")
dev.off()

#Print out filtered samples to file
ERBB2_RS_filter_samples=names(which(x>=ERBB2_RS_cutoff))
write.table(ERBB2_RS_filter_samples, file=outfile6B, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#Create new column to store ERBB2 status determined from array
ERBB2_RS_filtered=which(x >= ERBB2_RS_cutoff)
clin_data[,"ERBB2amp_array"]=0
clin_data[ERBB2_RS_filtered,"ERBB2amp_array"]=1

#Attempt to separate ESR1 data based on gaussian mixture model
#Allow mclust to determine best model, then filter out extreme dist'n (Want to filter ER-, therefore -> left/first component)
x=as.numeric(ESR1_Data[,4:length(header)])
mclust_ESR1=Mclust(x)
summary(mclust_ESR1, x) #Gives you list of values returned
classification_ESR1=mclust_ESR1$classification
num_clusters=mclust_ESR1$G
pdf(file=outfile5C)
par(mfrow=c(2,1), oma=c(2,2,2,2))
hist(x[classification_ESR1>1], xlim=c(2,15), col="blue", xlab="Log2 GCRMA value", main="ESR1 components to retain")
hist(x[classification_ESR1==1], xlim=c(2,15), col="red", xlab="Log2 GCRMA value", main="ESR1 components to filter")
title(main="Separation of ESR1- by model-based clustering", outer=TRUE)
ESR1_cutoff=max(x[classification_ESR1==1]) #Use minimum value of last cluster as cutoff point
ESR1_filtered_count=length(which(x <= ESR1_cutoff))
dev.off()


#Histogram - ESR1
pdf(file=outfile5D)
hist(as.numeric(ESR1_Data[,4:length(header)]), xlim=c(2,15), breaks=40, col="blue", main=paste("ESR1 (",ESR1_probe,") histogram for all 998 samples", sep=""), xlab="Log2 GCRMA value", ylab="Frequency")
abline(v=ESR1_cutoff, col="red")
legend("left", legend=c(paste("cutoff=",ESR1_cutoff,sep=""),paste("N=",ESR1_filtered_count,sep="")),bty="n")
dev.off()

#How many of the ESR1 filtered are actually from GSE11121 that reported 44/200 ER-?
#Answer 31/77 ESR1 filtered were from GSE11121 (aka Schmidt_2008)
ESR1_filtered=which(x <= ESR1_cutoff)
table(clin_data[ESR1_filtered,3])

#Print out filtered samples to file
x=as.matrix(ESR1_Data[,4:length(header)])
samples=colnames(x)
ESR1_filter_samples=samples[which(x <= ESR1_cutoff)]
write.table(ESR1_filter_samples, file=outfile6C, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#Create new column to store ER status determined from array
ESR1_filtered=which(x <= ESR1_cutoff)
clin_data[,"ER_array"]=1
clin_data[ESR1_filtered,"ER_array"]=0

#Write new file with ESR1 and ERBB2 status
write.table(clin_data, file=newclinfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




