#Load the appropriate libraries
library(affy)
library(gcrma)
library(genefilter)
library(randomForest)
library(ROCR)
require(Hmisc)
library("gplots") #For advanced heatmaps

#Set working directory for data files
setwd("/home/obig/Projects/Lymphoma_drugs/clustering_classification/CEL_FILES/")
varimp_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.varImps.pdf"
MDS_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.MDS.pdf"
heatmap_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.heatmap.pdf"
heatmap_topN_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.heatmap_topN.pdf"
case_pred_outfile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.CasePredictions.txt"
heatmap_FLgrade_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.heatmap_FLgrades.pdf"
heatmap_topN_FLgrade_pdffile="Lymphoma_GCBcells_vs_DLBCL_vs_FL.heatmap_topN_FLgrades.pdf"

Data=ReadAffy()

#Run GCRMA on Data to background correct, normalize and summarize expression data
gcrmaData=gcrma(Data)

#Set working data for results files
setwd("/home/obig/Projects/Lymphoma_drugs/clustering_classification")

#Write gcrma normalized data to file
#First, reduce number of decimal places
gcrma_values_formatted=format(exprs(gcrmaData), digits=5)
write.table(gcrma_values_formatted, file = "DLBCL_vs_FL_vs_normGCBcell_gcrma.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")

#Get sample names and expression data from gcrma expression object
pheno=pData(gcrmaData)
X=exprs(gcrmaData)

#Preliminary gene filtering might be a good idea. 
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation):
#At least 20% of samples should have raw intensity greater than 100
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_gcrmaData=gcrmaData[filt,]

#Create a vector of probe names after filtering for output later
var_names_filt=rownames(exprs(filt_gcrmaData))

#Create a vector for the class labels ('GC B-cell'=0 vs 'FL'=1 vs 'DLBCL'=2).
class_labels=vector(length=29)
class_labels[1:24]="FL"
class_labels[25:29]="GC B-cell"
class_labels[30:50]="DLBCL"
target=as.factor(class_labels)

#Get potential predictor variables
predictor_data=t(exprs(filt_gcrmaData))

#Run RandomForests
#rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 1000, proximity=TRUE, sampsize=c(24,5))
rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 10000, proximity=TRUE)

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
#MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"))
MDSplot(rf_output, target, k=2, main="RandomForests MDS plot of FL (F) vs DLBCL (D) vs GC B-cell (G)", xlab="", ylab="", pch=target_labels, palette=c("darkgreen", "yellowgreen", "blue"))
dev.off()

#Determine which cases were misclassified based on OOB testing and by how much
#Margins represent the fraction of votes for the correct class minus the fraction for the incorrect class
#If class is determined by majority rules, then cases with positive margins are correct and vice versa
margins=margin(rf_output, target)
#combine the patient ids with target/known class, predicted class, votes, and margins
case_predictions=cbind(as.vector(target),as.vector(rf_output$predicted),rf_output$votes,as.vector(margins))
colnames(case_predictions)=c("Target","Predicted","DLBCL_votes","FL_votes","GC_votes","Margins")
misclass_cases=case_predictions[case_predictions[,"Target"]!=case_predictions[,"Predicted"],]

#Write case predictions to file
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Create a heatmap for all data
patientcolors=class_labels
patientcolors[patientcolors=="FL"]="#9ACD32" #yellowgreen
patientcolors[patientcolors=="DLBCL"]="#006400" #darkgreen"
patientcolors[patientcolors=="GC B-cell"]="#0000FF" #blue

FL_grades=vector(length=24)
FL_grades_short=vector(length=24)
FL_grades=c("FLN0191_grade_3A","FLN0211_grade_1","FLN0231_grade_1","FLN0251_grade_3A","FLN0271_grade_1","FLN0291_grade_1","FLN0311_grade_1","FLN0331_grade_1","FLN0351_grade_1","FLN0371_grade_2","FLN0391_grade_1","FLN0411_grade_2","FLN0431_grade_2","FLN0454_grade_1","FLN0481_grade_1","FLN0511_grade_2","FLN0541_grade_2","FLN0561_grade_1","FLN0581_grade_3B","FLN0601_grade_2","FLN0621_grade_1","FLN0641_grade_1","FLN0681_grade_1","FLN0743_grade_1")
FL_grades_short=c("FL_grade_3A","FL_grade_1","FL_grade_1","FL_grade_3A","FL_grade_1","FL_grade_1","FL_grade_1","FL_grade_1","FL_grade_1","FL_grade_2","FL_grade_1","FL_grade_2","FL_grade_2","FL_grade_1","FL_grade_1","FL_grade_2","FL_grade_2","FL_grade_1","FL_grade_3B","FL_grade_2","FL_grade_1","FL_grade_1","FL_grade_1","FL_grade_1")
row_labels_FLgrades_short=c(FL_grades_short,class_labels[25:50])

x=as.matrix(predictor_data)
#rownames(x)=class_labels
pdf(file=heatmap_pdffile)
#heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=patientcolors, cexRow=0.8, cexCol=0.60)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main="Hierarchical clustering using all probes", labRow=class_labels, RowSideColors=patientcolors, labCol=FALSE, cexRow=0.8, cexCol=0.60, cex.main=0.8)
dev.off()

pdf(file=heatmap_FLgrade_pdffile)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main="Hierarchical clustering using top 30 predictors", labRow=row_labels_FLgrades_short, RowSideColors=patientcolors, labCol=FALSE, cexRow=0.8, cexCol=0.60)
dev.off()



#Create a heatmap for just the most discriminatory probes
topN=30
sorted_rf_importances=sort(rf_importances[,4])
varImp_topN=sorted_rf_importances[length(sorted_rf_importances)-(topN-1)]
topN_rf_importances=rf_importances[,4]>=varImp_topN
x=x[,topN_rf_importances]

pdf(file=heatmap_topN_pdffile)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main="Hierarchical clustering using top 30 predictors", labRow=class_labels, RowSideColors=patientcolors, labCol=colnames(x), cexRow=0.8, cexCol=0.60)
dev.off()

pdf(file=heatmap_topN_FLgrade_pdffile)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main="Hierarchical clustering using top 30 predictors", labRow=row_labels_FLgrades_short, RowSideColors=patientcolors, labCol=colnames(x), cexRow=0.8, cexCol=0.60)
dev.off()
