library(randomForest)
library(mclust)
require(Hmisc)
library(ROCR)
library("gplots")

#Set drug of interest
#drug="GSK_Tykerb" #Lapatinib
drug="BIBW2992"
#drug="AKT1.2.inhibitor"

outdir=paste("C:/Users/Obi/Documents/My Dropbox/drug_predictors/RFpredictors/singleDrugs/",drug,"/",sep="")
dir.create(outdir)
setwd(outdir)
outfile=paste("RFoutput.txt_",drug,".txt",sep="")
case_pred_outfile="CasePredictions.txt"
varimp_pdffile="varImps.pdf"
MDS_pdffile="MDS.pdf"
ROC_pdffile="ROC.pdf"
vote_dist_pdffile="vote_dist.pdf"
waterfallfile="waterfall.pdf"
heatmapfile="heatmap.pdf"
heatmapfile_scaled="heatmap_scaled.pdf"

#Import cell line data
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]


#Import combined predictor data
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/combined/BCCL_combined_data.txt"
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:2]
raw_data=raw_data_import[,3:length(colnames(raw_data_import))]

#Fix column names
colnames(raw_data)[which(colnames(raw_data)=="X184A1")]="184A1"
colnames(raw_data)[which(colnames(raw_data)=="X184B5")]="184B5"
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"

#Filter down to core set of cell lines - must have drug data and at least one other molecular profiling data type
core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")
data=raw_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#Set row names
rownames(data)=paste(raw_feat_data[,"DataType"],raw_feat_data[,"ID"], sep="__")

#drug response data
drugdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_10Feb.csv"
raw_drugdata_import=read.csv(drugdatafile)

drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Fix row names
rownames(drug_data)[which(rownames(drug_data)=="Hs578T")]="HS578T"

#Retrieve data for only libraries in core cell line set
drug_data_filt=drug_data[core_cell_lines,]

#Transform data to -log values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#For each drug of interest, find predictors associated with response and build predictor for drug response
#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]

mean_cutoff=mean(drug_data_interest, na.rm=TRUE)
drug_data_interest_NA=which(is.na(drug_data_interest))
resistants=which(drug_data_interest<=mean_cutoff)
sensitives=which(drug_data_interest>mean_cutoff)

response_class=vector(length=length(drug_data_interest))
response_class[drug_data_interest_NA]=NA
response_class[sensitives]="sensitive"
response_class[resistants]="resistant"

#Exclude libs where response_class=NA
nonNA=which(!is.na(response_class))
data_nonNA=data[,nonNA]
response_class_nonNA=response_class[nonNA]
cell_line_data_nonNA=cell_line_data[nonNA,]

#Make sure there are at least 10 libs classified as sensitive and resistant
num_sensitive=length(which(response_class_nonNA=="sensitive"))
num_resistant=length(which(response_class_nonNA=="resistant"))

target=as.factor(response_class_nonNA)


#########################################################################
#
#
#Need to fix handling of categorical variables before imputation. See RFimpute script
#
#
#########################################################################


#Impute missing values. So far, only cell lines where drug data is NA have been excluded. Predictor data still contains NAs
data_imputed=rfImpute(x=t(data_nonNA), y=target, ntree=300, iter=5)#Increase numbers of ntree and iter for final result?

#Down-sampling?
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

#Run RF on imputed data. Recall, you must extract all but first column where rfImpute inserts target there
#rf_output=randomForest(x=data_imputed[,2:length(colnames(data_imputed))], y=target, importance = TRUE, ntree = 10001, proximity=TRUE)
rf_output=randomForest(x=data_imputed[,2:length(colnames(data_imputed))], y=target, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=1-overall_error
class1_error=paste(rownames(confusion)[1],"error rate =",format(confusion[1,3]*100, digits=4), sep=" ")
class2_error=paste(rownames(confusion)[2],"error rate =",format(confusion[2,3]*100, digits=4), sep=" ")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity =",format(sensitivity, digits=4), sep=" ")
spec_out=paste("specificity =",format(specificity, digits=4), sep=" ")
err_out=paste("overall error rate =",format(overall_error, digits=4),sep=" ")
acc_out=paste("overall accuracy =",format(overall_accuracy, digits=4),sep=" ")
misclass_1=paste(confusion[1,2], "/", num_resistant, rownames(confusion)[1],"misclassed", sep=" ")
misclass_2=paste(confusion[2,1], "/", num_sensitive, rownames(confusion)[2],"misclassed", sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Var. Imp. (Gini) for top 30 predictors")
dev.off()

#Produce back-to-back histogram of vote distributions for Sensitive and Resistant
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"sensitive"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for cell lines classified by RF', axes=TRUE, ylab="Fraction votes (sensitive)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()

#Create ROC curve plot and calculate AUC
#Can use Sensitive vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling sensitive vs resistant
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
AUC_out=paste("AUC =",format(AUC, digits=5, scientific=FALSE), sep=" ")
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()

#Create side by side plot with MDS plot and performance stats
#Produce MDS plot
pdf(file=MDS_pdffile)
par(mfrow=c(1,2), mar=c(0.5,0.5,2,0), oma=c(4,4,4,4))
target_labels=as.vector(target)
target_labels[target_labels=="sensitive"]="S"
target_labels[target_labels=="resistant"]="R"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")

#Add stats
plot.new()
stats_legend=c(
"Classification Stats",
paste("n =",length(target), "cell lines", sep=" "),
sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2,AUC_out
)
legend("left", legend=stats_legend, bty="n")
title(main=drug, outer=TRUE)

dev.off()



#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2,AUC_out), file=outfile, append=TRUE)


#Create heatmap of top X predictors
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

top30=sort(rf_importances[,4], decreasing=TRUE)[1:30]
top30_data=data_nonNA[names(top30),cell_lines_sorted]

#Manually add any predictors of interest
#GSK_Tykerb
#interest_list=c("RPPA__HER2p1248","EA__ERBB2","EA__GRB7","WB__GRB7","RPPA__HER2","RS__EJ1132172__E3a_E4a__ERBB2","RS__ER77698__ER13a__GRB7","EA__PERLD1","U133A__ERBB2","EA__PPP1R1B","RS__EJ1132638__E3a_E4a__C17orf37","EA__C17orf37","U133A__GRB7","RS__G8360__GRB7__GRB7","RS__ER77610__ER6a__STARD3","RS__G8359__ERBB2__ERBB2","EA__PRSS8","EA__FLVCR2","WB__ERBB2","U133A__PGAP3","EA__CADPS2","CNV__ERBB2","RS__EJ1894321__E18a_E19a__FMNL1","EA__ST6GAL1","RS__G6636__STARD3__STARD3","RS__G8361__C17orf37__C17orf37","EA__BCAS1","EA__TGFB1",
#"RS__EB180791__E10_Ab__KHDC1","RS__ER77722__ER4b__PERLD1","RS__G6638__PPP1R1B__PPP1R1B","EA__7A5","EA__SPRY1","EA__GRAMD2","RS__ER113302__ER4a__CYP4X1","RS__ER140525__ER7a__NFATC2","EA__RAPGEF5","RS__G12505__PSCA__PSCA","U133A__ELF3","EA__ACY3","EA__TUBA1A","RS__EJ970734__E12b_E13a__SPOCD1","RS__ER110114__ER13a__SPOCD1","RS__ER246632__ER3a__PSCA","EA__ENO2","RS__ER77712__ER4a__PNMT","Meth__cg13263114__ERBB2","RS__EJ666293__E4a_E5a__QPCT","RS__ER77648__ER4a__PPP1R1B","EA__CEACAM6")

#AKT1.2.inhibitor
#interest_list=c("RS__EJ958319__E11a_E12a__IRAK2","RS__ER98367__ER12a__AXL","RS__EJ1229369__E18a_E19a__DENND2A","RS__EJ1800892__E4a_E5a__EFCAB4A","RS__ER260647__ER1a__LRRC26","Meth__cg02067021__DNAJC5B","RS__ER109771__ER14a__COL16A1","RS__EJ1888601__E11a_E12a__FAM70B","RS__ER26489__ER3a__AC009656.11","EA__GLS","U133A__CD44","RS__EJ1904990__E7a_E8a__ROR1","RS__ER8317__ER22a__ADAM12","RS__T41059__ENST00000375031__SH2D5","U133A__MT1E","RS__EJ1650854__E15a_E16a__C17orf28","EA__MSN","RS__ER177411__ER29a__CCDC88A","EA__CD44","EA__SLC16A1","RS__ER170150__ER22b__NRP2","EA__C17orf28","RS__ER274263__ER4a__SRPX","EA__AGPS","RS__ER194409__ER3a__NKX6-1","RS__EB352783__E6_Ab__CSF1",
#"Meth__cg10878307__ATAD4","RS__EJ1905118__E1a_E2b__MUC1","RS__ER48134__ER22a__LTBP2","EA__CREB3L4","RS__ER196490__ER6a__PRSS12","RS__EJ218462__E20a_E21a__SPEG","RS__ER83247__ER2a__SLC9A3R1","RS__ER105376__ER8a__C3","RS__EJ885219__E5a_E6a__C1orf38","EA__SPRY1","U133A__PTRF","U133A__SPDEF","RS__ER42761__ER18a__COL4A1","RS__ER267942__ER1a__CLCN4","RS__ER70311__ER4a__IL32","RS__EJ1116237__E2b_E3a__CMTM3","RS__EJ240034__E11a_E12a__ADD2","RS__ER137591__ER8a__SH2D5","RS__G21732__AL136362.10__AL136362.10","EA__SELENBP1","EA__C1orf64","RS__EB174773__E28_Aa__MICAL2","WB__SPDEF","RS__ER147030__ER15a__C21orf7")

#BIBW2992
interest_list=c("RS__ER122669__ER4a__ANP32E","RS__ER168389__ER18a__HECW2","RS__EJ913328__E1a_E2a__STARD3","RS__EJ1239593__E15a_E16a__WHSC1L1","ES__NSUN5C","RS__ER89231__ER5a__AC107983.3","ES__UGT2B11","RS__EJ1097536__E10a_E11a__ESD","RS__ER176017__ER3a__QPCT","RS__EJ171958__E16a_E17a__ANKRD44","RS__ER130538__ER12a__KIF21B","RS__EJ2094041__E4a_E5a__INPP5F","ES__SLC26A8","RS__ER151036__ER2a__SEPT3","RS__ER91904__ER6a__SLC39A6",
"RS__ER79547__ER5a__FMNL1","RS__EJ2052842__E9a_E10a__PLCG2","EA__RNF38","RS__EJ1033886__E2a_E3a__TUBB2B","RS__EJ1582024__E13a_E14a__DLC1","RS__EJ1132636__E2a_E3a__C17orf37","RS__G18828__AC008629.7__AC008629.7","RS__EJ1840763__E3a_E4a__AP1S2","RS__EJ1224112__E10a_E11a__SDK1","RS__ER53722__ER10a__SPG3A","RS__EJ2056198__E56a_E57a__VPS13A","RS__ER216234__ER2a__C6orf126","RS__ER10231__ER7a__GPR158","RS__ER72792__ER6a__SCNN1G","RS__EJ164095__E6a_E7a__HIPK2","WB__GRB7")


#top30_data=rbind(top30_data,data_nonNA["MM__PIK3CAmut",])
top30_data=data_nonNA[interest_list,cell_lines_sorted]

predictor_names=rownames(top30_data)
lib_names=colnames(top30_data)

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[cell_lines_sorted,"BCCLclassification2"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="ERBB2Amp"]="blue"
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"
subtype_colors[subtype_colors=="Non-malignant"]="yellow"
subtype_colors[subtype_colors=="Unknown"]="pink"

#Get data_type from predictor names
getPrefix=function(x){y=x[[1]][1]; return(y)}
predictor_names_split=strsplit(x=predictor_names, split="__")
prefixes=sapply(predictor_names_split,getPrefix, simplify = TRUE)

#Make predictor names without prefix, drop leading ID for RNAseq Names (only type to have 3 components)
getName=function(x){y=x[2:length(x)]; return(y)}
Names=sapply(predictor_names_split,getName)
pasteNames=function(x){
 if (length(x)==3){
  y=paste(x[2:3],collapse="__");
 }else{
  y=paste(x,collapse="__");
 }
return(y)
}
CleanNames=sapply(Names,pasteNames)

#Set colors for datatype sidebar
colors=rainbow(9)
datatype_colors=prefixes
datatype_colors[datatype_colors=="RS"]=colors[1]
datatype_colors[datatype_colors=="EA"]=colors[2]
datatype_colors[datatype_colors=="U133A"]=colors[3]
datatype_colors[datatype_colors=="ES"]=colors[4]
datatype_colors[datatype_colors=="MM"]=colors[5]
datatype_colors[datatype_colors=="CNV"]=colors[6]
datatype_colors[datatype_colors=="RPPA"]=colors[7]
datatype_colors[datatype_colors=="WB"]=colors[8]
datatype_colors[datatype_colors=="Meth"]=colors[9]

main_title=drug
pdf(file=heatmapfile)
heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=CleanNames, colsep=num_sensitive, sepcolor="black", cexRow=0.8, cexCol=0.75, margins=c(6,12), col=rev(heat.colors(75)))
legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","ExomeSeq","ManMut","CNV","RPPA","WB","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.7)
legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)
dev.off()

pdf(file=heatmapfile_scaled)
heatmap.2(as.matrix(top30_data), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="none", main=main_title, ColSideColors=subtype_colors, RowSideColors=datatype_colors, labCol=lib_names, labRow=CleanNames, colsep=num_sensitive, sepcolor="black", cexRow=0.7, cexCol=0.9, margins=c(6,12), col=rev(heat.colors(75)))
legend("bottomleft", legend=c("RNAseq","ExonArray","U133A","ExomeSeq","ManMut","CNV","RPPA","WB","Meth"), title="Data type", fill=colors, border="white", bty="n", cex=0.7)
legend("top", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), title="Subtype", border="white",  bty="n",cex=0.7)
dev.off()

#Create corresponding histogram/waterfall plot for drug response data
ymin=floor(min(drug_data_interest_sorted))
ymax=ceiling(max(drug_data_interest_sorted))
barplot_values=barplot(drug_data_interest_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
pdf(file=waterfallfile)
barplot(drug_data_interest_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
#draw line midway between bars straddling sensitive/resistant threshold
abline(v=barplot_values[num_sensitive]+(barplot_values[num_sensitive+1]-barplot_values[num_sensitive])/2, lwd=2, col="black")
legend("topright", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), cex=0.9)
dev.off()

