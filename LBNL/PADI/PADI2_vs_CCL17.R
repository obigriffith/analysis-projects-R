working_dir="/Users/ogriffit/Dropbox/LBNL/Projects/PADI/v3"
setwd(working_dir)


#Cell line info
celllinefile="/Users/ogriffit/Dropbox/LBNL/Projects/PADI/cell_line_info.2.txt"

#All gene expression data
datafile="/Users/ogriffit/Dropbox/LBNL/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="/Users/ogriffit/Dropbox/LBNL/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_cell_data=read.table(celllinefile, header = TRUE, na.strings = "NA", sep="\t")
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
#raw_gene_data=read.table(genedatafile, header = TRUE, na.strings = "NA", sep="\t")

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
high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"SubtypeNew"]!="Unknown") & as.vector(raw_cell_data[,"SubtypeNew"]!="Normal"))

#Apply library filter to datasets
data=raw_data[,high_qual]
feat_data=raw_feat_data
exp_status=raw_exp_status[,high_qual]
cell_data=raw_cell_data[high_qual,]
libs=colnames(data)

#Retrieve cell line details
lib_names=as.vector(cell_data[,"Sample.Name"])
subtypes=as.vector(cell_data[,"SubtypeNew"])
ERBB2=as.vector(cell_data[,"ERBB2New"])
Qualities=as.vector(cell_data[,"Quality"])
cbind(lib_names,libs)

#Sort cell lines by subtype then, name
lib_names_sorted=lib_names[order(cell_data[,"SubtypeNew"], cell_data[,"Sample.Name"])]
subtypes_sorted=subtypes[order(subtypes)]

#Get PADI2/CCL17 gene data, and sort according to subtype order determined above
PADI2_data_sorted=data[which(feat_data[,"Seq_Name"]=="PADI2"),lib_names_sorted]
CCL17_data_sorted=data[which(feat_data[,"Seq_Name"]=="CCL17"),lib_names_sorted]

#Also get expression status for these genes and summarize
PADI2_exp_status_sorted=exp_status[which(feat_data[,"Seq_Name"]=="PADI2"),lib_names_sorted]
CCL17_exp_status_sorted=exp_status[which(feat_data[,"Seq_Name"]=="CCL17"),lib_names_sorted]

paste(sum(PADI2_exp_status_sorted[subtypes_sorted=="Basal"]), "/", length(which(subtypes_sorted=="Basal")), "Basal Expressed")
paste(sum(PADI2_exp_status_sorted[subtypes_sorted=="Basal_NM"]), "/", length(which(subtypes_sorted=="Basal_NM")), "Basal_NM Expressed")
paste(sum(PADI2_exp_status_sorted[subtypes_sorted=="ClaudinLow"]), "/", length(which(subtypes_sorted=="ClaudinLow")), "ClaudinLow Expressed")
paste(sum(PADI2_exp_status_sorted[subtypes_sorted=="Luminal"]), "/", length(which(subtypes_sorted=="Luminal")), "Luminal Expressed")

paste(sum(CCL17_exp_status_sorted[subtypes_sorted=="Basal"]), "/", length(which(subtypes_sorted=="Basal")), "Basal Expressed")
paste(sum(CCL17_exp_status_sorted[subtypes_sorted=="Basal_NM"]), "/", length(which(subtypes_sorted=="Basal_NM")), "Basal_NM Expressed")
paste(sum(CCL17_exp_status_sorted[subtypes_sorted=="ClaudinLow"]), "/", length(which(subtypes_sorted=="ClaudinLow")), "ClaudinLow Expressed")
paste(sum(CCL17_exp_status_sorted[subtypes_sorted=="Luminal"]), "/", length(which(subtypes_sorted=="Luminal")), "Luminal Expressed")



#Convert values to log2, unless they are 0 in which case set them to 0
PADI2_data_log2 = log2(PADI2_data_sorted+1)
CCL17_data_log2 = log2(CCL17_data_sorted+1)

#Calculate correlations between PADI2 and CCL17
cors=matrix(NA, nrow=4, ncol=2, dimnames = list(unique(subtypes_sorted),c("Spearman_rho", "Spearman_p")))
spearman_all=cor.test(x=as.numeric(CCL17_data_log2), y=as.numeric(PADI2_data_log2), method="spearman")
spearman_Basal=cor.test(x=as.numeric(CCL17_data_log2[which(subtypes_sorted=="Basal")]), y=as.numeric(CCL17_data_log2[which(subtypes_sorted=="Basal")]), method="spearman")
spearman_Basal_NM=cor.test(x=as.numeric(CCL17_data_log2[which(subtypes_sorted=="Basal_NM")]), y=as.numeric(PADI2_data_log2[which(subtypes_sorted=="Basal_NM")]), method="spearman")
spearman_Luminal=cor.test(x=as.numeric(CCL17_data_log2[which(subtypes_sorted=="Luminal")]), y=as.numeric(PADI2_data_log2[which(subtypes_sorted=="Luminal")]), method="spearman")
spearman_Claudin=cor.test(x=as.numeric(CCL17_data_log2[which(subtypes_sorted=="ClaudinLow")]), y=as.numeric(PADI2_data_log2[which(subtypes_sorted=="ClaudinLow")]), method="spearman")
spearman_notBasal=cor.test(x=as.numeric(CCL17_data_log2[which(subtypes_sorted!="Basal")]), y=as.numeric(PADI2_data_log2[which(subtypes_sorted!="Basal")]), method="spearman")


cors["Basal",]=c(spearman_Basal$estimate,spearman_Basal$p.value)
cors["Basal_NM",]=c(spearman_Basal_NM$estimate,spearman_Basal_NM$p.value)
cors["ClaudinLow",]=c(spearman_Claudin$estimate,spearman_Claudin$p.value)
cors["Luminal",]=c(spearman_Luminal$estimate,spearman_Luminal$p.value)

pdf(file="PADI2_CCL17_cor_by_subtype.pdf", width=10, height=7.5)
#Break plot into uneven sections for each subtype to be plotted proportionally
layout(matrix(c(1,1,1,1,1,1,2,2,3,3,4,4,4,4,4,4,4,4), 1, 18, byrow = TRUE))
#layout.show(4)
par(mar=c(7,4,3,1)) #bottom, left, top, right

#Create plot to represent correlation between PADI2 and CCL17
#break into panels for each subtype
i=1
subtype=unique(subtypes_sorted)[i]
CCL17_data_log2_sub=CCL17_data_log2[subtypes_sorted==subtype]
PADI2_data_log2_sub=PADI2_data_log2[subtypes_sorted==subtype]
lib_names_sub=names(CCL17_data_log2_sub)
plot(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), type="n", xlim=c(1,length(lib_names_sub)), ylim=c(0,max(c(as.numeric(CCL17_data_log2),as.numeric(PADI2_data_log2)))+1), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", bty="n")
axis(1, at=1:length(lib_names_sub), labels=lib_names_sub, las=2, cex.axis=0.9)
points(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), col="blue", pch=18, cex=1.6)
points(x=(1:length(lib_names_sub)), y=as.numeric(PADI2_data_log2_sub), col="red", pch=15, cex=1.2)
lines(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), col="blue", lwd=2)
lines(x=(1:length(lib_names_sub)), y=as.numeric(PADI2_data_log2_sub), col="red", lwd=2)
legend("top", legend=c(subtype, paste("rho =",format(cors[subtype,"Spearman_rho"], digit=3)),paste("p =",format(cors[subtype,"Spearman_p"],digit=3))), bty="n")

par(mar=c(7,2,3,1)) #bottom, left, top, right
for (i in 2:length(unique(subtypes_sorted))){
 subtype=unique(subtypes_sorted)[i]
 CCL17_data_log2_sub=CCL17_data_log2[subtypes_sorted==subtype]
 PADI2_data_log2_sub=PADI2_data_log2[subtypes_sorted==subtype]
 lib_names_sub=names(CCL17_data_log2_sub)
 plot(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), type="n", xlim=c(1,length(lib_names_sub)), ylim=c(0,max(c(as.numeric(CCL17_data_log2),as.numeric(PADI2_data_log2)))+1), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", yaxt="n", bty="n")
 axis(1, at=1:length(lib_names_sub), labels=lib_names_sub, las=2, cex.axis=0.9)
 points(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), col="blue", pch=18, cex=1.6)
 points(x=(1:length(lib_names_sub)), y=as.numeric(PADI2_data_log2_sub), col="red", pch=15, cex=1.2)
 lines(x=(1:length(lib_names_sub)), y=as.numeric(CCL17_data_log2_sub), col="blue", lwd=2)
 lines(x=(1:length(lib_names_sub)), y=as.numeric(PADI2_data_log2_sub), col="red", lwd=2)
 legend("top", legend=c(subtype, paste("rho =",format(cors[subtype,"Spearman_rho"], digit=3)),paste("p =",format(cors[subtype,"Spearman_p"],digit=3))), bty="n")
}
legend("bottomright", legend=c("CCL17","PADI2"), pch=c(18,15), lwd=c(2,2), col=c("blue","red"), pt.cex=c(1.5,1.2), bty="n")
dev.off()

####Do same thing for HER2+ versus HER2-#####
#Sort cell lines by ERBB2 then, then name
lib_names_sorted2=lib_names[order(cell_data[,"ERBB2New"], cell_data[,"Sample.Name"])]
ERBB2_sorted=ERBB2[order(ERBB2)]

#Get PADI2/CCL17 gene data, and sort according to ERBB2 order determined above
CCL17_data_sorted2=data[which(feat_data[,"Seq_Name"]=="CCL17"),lib_names_sorted2]
PADI2_data_sorted2=data[which(feat_data[,"Seq_Name"]=="PADI2"),lib_names_sorted2]

#Convert values to log2, unless they are 0 in which case set them to 0
CCL17_data_sorted2_log2 = log2(CCL17_data_sorted2+1)
PADI2_data_sorted2_log2 = log2(PADI2_data_sorted2+1)

#Calculate correlations between PADI2 and CCL17
spearman_ERBB2pos=cor.test(x=as.numeric(CCL17_data_sorted2_log2[which(ERBB2_sorted=="Amp")]), y=as.numeric(PADI2_data_sorted2_log2[which(ERBB2_sorted=="Amp")]), method="spearman")
spearman_ERBB2neg=cor.test(x=as.numeric(CCL17_data_sorted2_log2[which(ERBB2_sorted=="NoAmp")]), y=as.numeric(PADI2_data_sorted2_log2[which(ERBB2_sorted=="NoAmp")]), method="spearman")

pdf(file="PADI2_CCL17_cor_by_ERBB2.pdf", width=10, height=7.5)
#Break plot into sections for each subtype to be plotted proportionally
layout(matrix(c(1,1,1,1,2,2,2,2,2,2), 1, 10, byrow = TRUE))
#layout.show(2)

par(mar=c(7,4,3,1)) #bottom, left, top, right
#Create plot to represent correlation between PADI2 and CCL17
#break into panels for each ERBB2 category
CCL17_data_log2_ERBB2pos=CCL17_data_sorted2_log2[ERBB2_sorted=="Amp"]
PADI2_data_log2_ERBB2pos=PADI2_data_sorted2_log2[ERBB2_sorted=="Amp"]
lib_names_ERBB2pos=names(CCL17_data_log2_ERBB2pos)
plot(x=(1:length(lib_names_ERBB2pos)), y=as.numeric(CCL17_data_log2_ERBB2pos), type="n", xlim=c(1,length(lib_names_ERBB2pos)), ylim=c(0,max(c(as.numeric(CCL17_data_log2),as.numeric(PADI2_data_log2)))+1), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", bty="n")
axis(1, at=1:length(lib_names_ERBB2pos), labels=lib_names_ERBB2pos, las=2, cex.axis=0.9)
points(x=(1:length(lib_names_ERBB2pos)), y=as.numeric(CCL17_data_log2_ERBB2pos), col="blue", pch=18, cex=1.6)
points(x=(1:length(lib_names_ERBB2pos)), y=as.numeric(PADI2_data_log2_ERBB2pos), col="red", pch=15, cex=1.2)
lines(x=(1:length(lib_names_ERBB2pos)), y=as.numeric(CCL17_data_log2_ERBB2pos), col="blue", lwd=2)
lines(x=(1:length(lib_names_ERBB2pos)), y=as.numeric(PADI2_data_log2_ERBB2pos), col="red", lwd=2)
legend("top", legend=c("ERBB2+", paste("rho =",format(spearman_ERBB2pos$estimate, digit=3)),paste("p =",format(spearman_ERBB2pos$p.value,digit=3))), bty="n")

par(mar=c(7,2,3,1)) #bottom, left, top, right
CCL17_data_log2_ERBB2neg=CCL17_data_sorted2_log2[ERBB2_sorted=="NoAmp"]
PADI2_data_log2_ERBB2neg=PADI2_data_sorted2_log2[ERBB2_sorted=="NoAmp"]
lib_names_ERBB2neg=names(CCL17_data_log2_ERBB2neg)
plot(x=(1:length(lib_names_ERBB2neg)), y=as.numeric(CCL17_data_log2_ERBB2neg), type="n", xlim=c(1,length(lib_names_ERBB2neg)), ylim=c(0,max(c(as.numeric(CCL17_data_log2),as.numeric(PADI2_data_log2)))+1), xlab=NA, ylab="Log2 Gene Expression Level", xaxt="n", yaxt="n", bty="n")
axis(1, at=1:length(lib_names_ERBB2neg), labels=lib_names_ERBB2neg, las=2, cex.axis=0.9)
points(x=(1:length(lib_names_ERBB2neg)), y=as.numeric(CCL17_data_log2_ERBB2neg), col="blue", pch=18, cex=1.6)
points(x=(1:length(lib_names_ERBB2neg)), y=as.numeric(PADI2_data_log2_ERBB2neg), col="red", pch=15, cex=1.2)
lines(x=(1:length(lib_names_ERBB2neg)), y=as.numeric(CCL17_data_log2_ERBB2neg), col="blue", lwd=2)
lines(x=(1:length(lib_names_ERBB2neg)), y=as.numeric(PADI2_data_log2_ERBB2neg), col="red", lwd=2)
legend("top", legend=c("ERBB2-", paste("rho =",format(spearman_ERBB2neg$estimate, digit=3)),paste("p =",format(spearman_ERBB2neg$p.value,digit=3))), bty="n")
legend("bottomright", legend=c("ERBB2","PADI2"), pch=c(18,15), lwd=c(2,2), col=c("blue","red"), pt.cex=c(1.5,1.2), bty="n")
dev.off()




#CODE BELOW TO BE UPDATED FOR CCL17 instead of ERBB2 if needed


###Do same thing for Basal/ERBB2+ vs Basal/ERBB2- vs Luminal/ERBB2+ vs Luminal ERBB2-###


#Produce box and whisker plot
#PADI2
pdf(file="PADI2_expression_all_lines.pdf", width=10, height=7.5)
subtype_rainbow=rainbow(length(unique(subtypes_sorted)))
names(subtype_rainbow)=unique(subtypes_sorted)
subtype_colors=subtype_rainbow[subtypes_sorted]

#Replace zero values with NA. Works better with boxplot (and is consistent with alexa-seq plots)
data_noZero=data[,lib_names_sorted]
data_noZero[data_noZero==0]=NA
PADI2_data_noZero=PADI2_data_sorted
PADI2_data_noZero[PADI2_data_noZero==0]=NA

#Histogram of gene expression data showing position of current gene
main_title = "PADI2 gene-level expression compared to distribution of all genes"
y_label = "log2 gene expression level"
boxplot(x = log2(data_noZero), col=subtype_colors, main=main_title, ylab=y_label, las=2, col.lab = gray(.1), cex.main = 1, cex.lab = 1, cex.axis=0.7)
target_gene_level = as.numeric(log2(PADI2_data_noZero))
points(target_gene_level, pch=16, col="black", cex=1)
legend("topleft", legend=c(unique(subtypes_sorted),"PADI2"), fill=c(unique(subtype_colors),"white"), pch=c(NA,NA,NA,NA,16), bg="white", border=NA)
dev.off()





#PADI4
pdf(file="PADI4_expression_all_lines.pdf", width=10, height=7.5)
subtype_rainbow=rainbow(length(unique(subtypes_sorted)))
names(subtype_rainbow)=unique(subtypes_sorted)
subtype_colors=subtype_rainbow[subtypes_sorted]

#Replace zero values with NA. Works better with boxplot (and is consistent with alexa-seq plots)
data_noZero=data[,lib_names_sorted]
data_noZero[data_noZero==0]=NA
PADI4_data_noZero=PADI4_data_sorted
PADI4_data_noZero[PADI4_data_noZero==0]=NA

#Histogram of gene expression data showing position of current gene
main_title = "PADI4 gene-level expression compared to distribution of all genes"
y_label = "log2 gene expression level"
boxplot(x = log2(data_noZero), col=subtype_colors, main=main_title, ylab=y_label, las=2, col.lab = gray(.1), cex.main = 1, cex.lab = 1, cex.axis=0.7)
target_gene_level = as.numeric(log2(PADI4_data_noZero))
points(target_gene_level, pch=16, col="black", cex=1)
legend("topleft", legend=c(unique(subtypes_sorted),"PADI4"), fill=c(unique(subtype_colors),"white"), pch=c(NA,NA,NA,NA,16), bg="white", border=NA)
dev.off()


#Calculate statistics for PADI2 differential expression between subtypes and ERBB2 status
###Subtype statistics###
PADI2_subtype_data=list(Basal=as.numeric(PADI2_data_log2[subtypes_sorted=="Basal"]), Basal_NM=as.numeric(PADI2_data_log2[subtypes_sorted=="Basal_NM"]), ClaudinLow=as.numeric(PADI2_data_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(PADI2_data_log2[subtypes_sorted=="Luminal"])) 
PADI4_subtype_data=list(Basal=as.numeric(PADI4_data_log2[subtypes_sorted=="Basal"]), Basal_NM=as.numeric(PADI4_data_log2[subtypes_sorted=="Basal_NM"]), ClaudinLow=as.numeric(PADI4_data_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(PADI4_data_log2[subtypes_sorted=="Luminal"])) 
ERBB2_subtype_data=list(Basal=as.numeric(ERBB2_data_log2[subtypes_sorted=="Basal"]), Basal_NM=as.numeric(ERBB2_data_log2[subtypes_sorted=="Basal_NM"]), ClaudinLow=as.numeric(ERBB2_data_log2[subtypes_sorted=="ClaudinLow"]), Luminal=as.numeric(ERBB2_data_log2[subtypes_sorted=="Luminal"])) 

#PADI2
#basal vs luminal
wilcox.test(x=PADI2_subtype_data$Basal, y=PADI2_subtype_data$Luminal, alternative="two.sided")
#basal vs claudinlow/basal_nm
wilcox.test(x=PADI2_subtype_data$Basal, y=c(PADI2_subtype_data$Basal_NM,PADI2_subtype_data$ClaudinLow), alternative="two.sided")
#basal/luminal vs claudinlow/basal_nm
wilcox.test(x=c(PADI2_subtype_data$Basal,PADI2_subtype_data$Luminal), y=c(PADI2_subtype_data$Basal_NM,PADI2_subtype_data$ClaudinLow), alternative="two.sided")
#basal_nm vs claudinlow
wilcox.test(x=PADI2_subtype_data$Basal_NM, y=PADI2_subtype_data$ClaudinLow, alternative="two.sided")
#luminal vs basal/basal_nm/claudin-low
wilcox.test(x=PADI2_subtype_data$Luminal, y=c(PADI2_subtype_data$Basal,PADI2_subtype_data$Basal_NM,PADI2_subtype_data$ClaudinLow), alternative="two.sided")

#PADI4
#basal vs luminal
wilcox.test(x=PADI4_subtype_data$Basal, y=PADI4_subtype_data$Luminal, alternative="two.sided")
#basal vs claudinlow/basal_nm
wilcox.test(x=PADI4_subtype_data$Basal, y=c(PADI4_subtype_data$Basal_NM,PADI4_subtype_data$ClaudinLow), alternative="two.sided")
#basal/luminal vs claudinlow/basal_nm
wilcox.test(x=c(PADI4_subtype_data$Basal,PADI4_subtype_data$Luminal), y=c(PADI4_subtype_data$Basal_NM,PADI4_subtype_data$ClaudinLow), alternative="two.sided")
#basal_nm vs claudinlow
wilcox.test(x=PADI4_subtype_data$Basal_NM, y=PADI4_subtype_data$ClaudinLow, alternative="two.sided")
#luminal vs basal/basal_nm/claudin-low
wilcox.test(x=PADI4_subtype_data$Luminal, y=c(PADI4_subtype_data$Basal,PADI4_subtype_data$Basal_NM,PADI4_subtype_data$ClaudinLow), alternative="two.sided")
#basal_nm vs basal/luminal/claudin-low
wilcox.test(x=PADI4_subtype_data$Basal_NM, y=c(PADI4_subtype_data$Basal,PADI4_subtype_data$Luminal,PADI4_subtype_data$ClaudinLow), alternative="two.sided")





#ERBB2
#basal vs luminal
wilcox.test(x=ERBB2_subtype_data$Basal, y=ERBB2_subtype_data$Luminal, alternative="two.sided")
#basal vs claudinlow/basal_nm
wilcox.test(x=ERBB2_subtype_data$Basal, y=c(ERBB2_subtype_data$Basal_NM,ERBB2_subtype_data$ClaudinLow), alternative="two.sided")
#basal/luminal vs claudinlow/basal_nm
wilcox.test(x=c(ERBB2_subtype_data$Basal,ERBB2_subtype_data$Luminal), y=c(ERBB2_subtype_data$Basal_NM,ERBB2_subtype_data$ClaudinLow), alternative="two.sided")
#basal_nm vs claudinlow
wilcox.test(x=ERBB2_subtype_data$Basal_NM, y=ERBB2_subtype_data$ClaudinLow, alternative="two.sided")

#Boxplot of PADI2 by subtype
pdf(file="PADI2_subtype_comp.pdf", width=7.5, height=7.5)
ymin=min(as.numeric(PADI2_data_log2))
ymax=max(as.numeric(PADI2_data_log2))
par(mar=c(8,5,3,1)) #bottom, left, top, right
boxplot(x=PADI2_subtype_data, las=2, col=subtype_rainbow, main="PADI2", ylab=y_label, ylim=c(ymin,ymax), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

#Boxplot of PADI4 by subtype
pdf(file="PADI4_subtype_comp.pdf", width=7.5, height=7.5)
ymin=min(as.numeric(PADI4_data_log2))
ymax=max(as.numeric(PADI4_data_log2))
par(mar=c(8,5,3,1)) #bottom, left, top, right
boxplot(x=PADI4_subtype_data, las=2, col=subtype_rainbow, main="PADI4", ylab=y_label, ylim=c(ymin,ymax), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()


#Boxplots of PADI2 and ERBB2 by subtype
pdf(file="PADI2_ERBB2_subtype_comp.pdf", width=7.5, height=7.5)
ymin=min(c(as.numeric(PADI2_data_log2),as.numeric(ERBB2_data_log2)))
ymax=max(c(as.numeric(PADI2_data_log2),as.numeric(ERBB2_data_log2)))
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar=c(7,4,3,1)) #bottom, left, top, right
boxplot(x=PADI2_subtype_data, las=2, col=subtype_rainbow, main="PADI2", ylab=y_label, ylim=c(ymin,ymax))
boxplot(x=ERBB2_subtype_data, las=2, col=subtype_rainbow, main="ERBB2", ylab=y_label, ylim=c(ymin,ymax))
dev.off()

###ERBB2 statistics###
PADI2_ERBB2_data=list(ERBB2pos=as.numeric(PADI2_data_sorted2_log2[ERBB2_sorted=="Amp"]), ERBB2neg=as.numeric(PADI2_data_sorted2_log2[ERBB2_sorted=="NoAmp"])) 
ERBB2_ERBB2_data=list(ERBB2pos=as.numeric(ERBB2_data_sorted2_log2[ERBB2_sorted=="Amp"]), ERBB2neg=as.numeric(ERBB2_data_sorted2_log2[ERBB2_sorted=="NoAmp"])) 

#PADI2
#ERBB2pos vs ERBB2neg
PADI2_ERBB2amp_stats=wilcox.test(x=PADI2_ERBB2_data$ERBB2pos, y=PADI2_ERBB2_data$ERBB2neg, alternative="two.sided")

#ERBB2
#ERBB2pos vs ERBB2neg
ERBB2_ERBB2amp_stats=wilcox.test(x=ERBB2_ERBB2_data$ERBB2pos, y=ERBB2_ERBB2_data$ERBB2neg, alternative="two.sided")

pdf(file="PADI2_ERBB2_ERBB2amp_comp.pdf", width=7.5, height=7.5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar=c(7,4,3,1)) #bottom, left, top, right
boxplot(x=PADI2_ERBB2_data, las=2, col=subtype_rainbow, main="PADI2", ylab=y_label, ylim=c(ymin,ymax))
boxplot(x=ERBB2_ERBB2_data, las=2, col=subtype_rainbow, main="ERBB2", ylab=y_label, ylim=c(ymin,ymax))
dev.off()


#Determine distribution of correlations between PADI2 and all genes
#Get PADI2 gene data, and sort according to subtype order determined above
PADI2_data_sorted=data[which(feat_data[,"Seq_Name"]=="PADI2"),lib_names_sorted]
nonPADI2_data_sorted=data[which(feat_data[,"Seq_Name"]!="PADI2"),lib_names_sorted]
gene_names=feat_data[,"Seq_Name"]

#Convert values to log2, unless they are 0 in which case set them to 0
PADI2_data_log2 = log2(PADI2_data_sorted+1)
nonPADI2_data_log2 = log2(nonPADI2_data_sorted+1)

cors=rep(NA, length(rownames(nonPADI2_data_log2)))
for (i in 1:length(rownames(nonPADI2_data_log2))){
 cor=cor.test(x=as.numeric(PADI2_data_log2), y=as.numeric(nonPADI2_data_log2[i,]), method="spearman")
 cors[i]=cor$estimate 
}
length(cors[which(cors>spearman_all$estimate)]) #Determine number of genes with better correlation than PADI2/ERBB2


#Get ERBB2 gene data, and sort according to subtype order determined above
ERBB2_data_sorted=data[which(feat_data[,"Seq_Name"]=="ERBB2"),lib_names_sorted]
nonERBB2_data_sorted=data[which(feat_data[,"Seq_Name"]!="ERBB2"),lib_names_sorted]

#Convert values to log2, unless they are 0 in which case set them to 0
ERBB2_data_log2 = log2(ERBB2_data_sorted+1)
nonERBB2_data_log2 = log2(nonERBB2_data_sorted+1)

cors=rep(NA, length(rownames(nonERBB2_data_log2)))
for (i in 1:length(rownames(nonERBB2_data_log2))){
 cor=cor.test(x=as.numeric(ERBB2_data_log2), y=as.numeric(nonERBB2_data_log2[i,]), method="spearman")
 cors[i]=cor$estimate 
}
length(cors[which(cors>spearman_all$estimate)]) #Determine number of genes with better correlation than PADI2/ERBB2


#Exclude basal from data, as these are actually negatively correlated in PADI2/ERBB2 correlation
#PADI2 vs all genes
PADI2_data_sorted_nonBasal=PADI2_data_sorted[subtypes_sorted!="Basal"]
nonPADI2_data_sorted_nonBasal=nonPADI2_data_sorted[,subtypes_sorted!="Basal"]

#Convert values to log2, unless they are 0 in which case set them to 0
PADI2_data_nonBasal_log2 = log2(PADI2_data_sorted_nonBasal+1)
nonPADI2_data_nonBasal_log2 = log2(nonPADI2_data_sorted_nonBasal+1)

cors2=rep(NA, length(rownames(nonPADI2_data_nonBasal_log2)))
for (i in 1:length(rownames(nonPADI2_data_nonBasal_log2))){
 cor=cor.test(x=as.numeric(PADI2_data_nonBasal_log2), y=as.numeric(nonPADI2_data_nonBasal_log2[i,]), method="spearman")
 cors2[i]=cor$estimate 
}
length(cors2[which(cors2>spearman_notBasal$estimate)])

#ERBB2 vs all genes
ERBB2_data_sorted_nonBasal=ERBB2_data_sorted[subtypes_sorted!="Basal"]
nonERBB2_data_sorted_nonBasal=nonERBB2_data_sorted[,subtypes_sorted!="Basal"]

#Convert values to log2, unless they are 0 in which case set them to 0
ERBB2_data_nonBasal_log2 = log2(ERBB2_data_sorted_nonBasal+1)
nonERBB2_data_nonBasal_log2 = log2(nonERBB2_data_sorted_nonBasal+1)

cors2=rep(NA, length(rownames(nonERBB2_data_nonBasal_log2)))
for (i in 1:length(rownames(nonERBB2_data_nonBasal_log2))){
 cor=cor.test(x=as.numeric(ERBB2_data_nonBasal_log2), y=as.numeric(nonERBB2_data_nonBasal_log2[i,]), method="spearman")
 cors2[i]=cor$estimate 
}
length(cors2[which(cors2>spearman_notBasal$estimate)]) #Number of genes which correlated with ERBB2 better than PADI2
ERBB2_better_cors=cors2[which(cors2>spearman_notBasal$estimate)] #correlation values of genes which correlated with ERBB2 better than PADI2
ERBB2_better_cors_names=feat_data[,"Seq_Name"][which(cors2>spearman_notBasal$estimate)] #names of genes which correlated with ERBB2 better than PADI2
cbind(ERBB2_better_cors_names,ERBB2_better_cors)


PADI2_data_sorted=data[which(feat_data[,"Seq_Name"]=="PADI2"),lib_names_sorted]
CCL17_data_sorted=data[which(feat_data[,"Seq_Name"]=="CCL17"),lib_names_sorted]
ERBB2_data_sorted=data[which(feat_data[,"Seq_Name"]=="ERBB2"),lib_names_sorted]



cor(x=as.numeric(PADI2_data_sorted[subtypes_sorted!="Basal"]),y=as.numeric(CCL17_data_sorted[subtypes_sorted!="Basal"]), method="spearman")


