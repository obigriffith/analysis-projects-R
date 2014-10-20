
#Load SNP6 segmented data
SNP6_segmented_data_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/CNV/prefiltered/JWGray_BCCL_SNP6_segmented_v2_table.txt"
SNP6_segmented_data=read.table(SNP6_segmented_data_file, header=TRUE)
SNP6_segments=SNP6_segmented_data[,1:3]
SNP6_segment_values=SNP6_segmented_data[,4:length(colnames(SNP6_segmented_data))]

#Fix colnames
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X184B5"]="184B5"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X185A1"]="185A1"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X600MPE"]="600MPE"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X21MT1"]="21MT1"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X21MT2"]="21MT2"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X21NT"]="21NT"
colnames(SNP6_segment_values)[colnames(SNP6_segment_values)=="X21PT"]="21PT"

#Load SNP6 GISTIC data
GISTIC_data_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/CNV/GISTIC_ROI_10JunV2.clean.2.txt"
GISTIC_data=read.table(GISTIC_data_file, header=TRUE)

#Pull out common set of cell lines (all GISTIC lines are in SNP6 segmented set)
SNP6_segment_values=SNP6_segment_values[,rownames(GISTIC_data)]

gene="PIK3CA"; chrom=3; gene_start=180349005; gene_end=180435191
gene="MYC"; chrom=8; gene_start=128817498; gene_end=128822855
gene="CCND1"; chrom=11; gene_start=69165054; gene_end=69178423
gene="ERBB2"; chrom=17; gene_start=35109780; gene_end=35138441
gene="AURKA"; chrom=20; gene_start=54377852; gene_end=54400758
gene="CDKN2A"; chrom=9; gene_start=21957751; gene_end=21984490
gene="CDKN2B"; chrom=9; gene_start=21992902; gene_end=21999312
gene="PTEN"; chrom=10; gene_start=89613175; gene_end=89718512

pdf(file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/CNV/SNP6_segmented_vs_GISTIC.pdf")
par(mfrow=c(3,3))

gene_segments=which(SNP6_segments[,"chrom"]==chrom & SNP6_segments[,"end"]>gene_start & SNP6_segments[,"start"]<gene_end)
gene_values=SNP6_segment_values[gene_segments,]
gene_values_mean=mean(gene_values, na.rm=TRUE)
#Creat plots to show difference in mean segmented values according to GISTIC groups
colors=rainbow(3)
boxplot(gene_values_mean[GISTIC_data[,"ERBB2"]==0], gene_values_mean[GISTIC_data[,"ERBB2"]==1], gene_values_mean[GISTIC_data[,"ERBB2"]==2], col=colors, names=c("0","1","2"), xlab="GISTIC", ylab="mean log2 rel. copy num. (seg)", main=gene)


dev.off()



