library("gplots") #For advanced heatmaps

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/processed_marker_data/62markers_23JAN07/clustering_and_classification")

data=read.table("ATC_all_deconvoluted_data_62markers_23JAN07_matched_diff_vs_undiff.txt", header = TRUE, na.strings = "NA", sep="\t")

marker_data=data[,c(6:67)]

pathology=data[,2]
pathology[pathology==1]="A"
pathology[pathology==3]="F"
pathology2=data[,5]
pathology[pathology2=="C-papillary"]="P"
pathology_patient=paste(pathology,data[,1], sep="")
patient_numbers_new=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
pathology_patient_new=paste(pathology,patient_numbers_new, sep="")


#Create a vector assigning a color to each patient based on pathology for the color side bar
patientcolors=pathology
patientcolors[patientcolors=="A"]="#9ACD32" #yellowgreen
#patientcolors[patientcolors=="F"]="#FF0000" #red
patientcolors[patientcolors=="F"]="#006400" #darkgreen
patientcolors[patientcolors=="P"]="#006400" #darkgreen

#To create a heatmap of all data
#Try different color scheme and add a score key
#pdf("ATC_all_deconvoluted_data_62markers_23JAN07_matched_heatmap_wkey.pdf")
pdf("ATC_all_deconvoluted_data_62markers_23JAN07_matched_heatmap_wkey.2.pdf")
#Choose colors for the five possible scores (0,1,2,3,4)
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names separately
#all_col_names=c("AAT","AMF","Aurora-A","Aurora-B","Aurora-C","BCL2","B-Catenin","CA9","CAV1","CDX2","CKIT","CR3","Cyclin-D1","Cyclin-E","Caveolin","Clusterin","ECAD","EGFR","ER","HER2","HER3","HER4","HSP27","IGF1R","ILK","INH","MDM2","MIB1","O13","P16","P21","P27","P53","P57","P63","P-AKT","PGI","PR","PSA","AR","CK19","COX2","Galectin","HBME-1","MLH1","MSH2","MSH6","P504S","P75-NTR","PMS2","TOPO2","TS106","RET","S100","Syntrophin","TDT","TSH","TTF1","Thyroglobulin","UPA","VEGF","WT1")
#all_col_names=c("AAT","AMF-R","Aurora-A","Aurora-B","Aurora-C","Bcl-2","CTNNB1","CA9","CAV-1","CDX2","c-kit","CR3","Cyclin-D1","Cyclin-E","Caveolin","Clusterin","E-CAD","EGFR","ER","HER2","HER3","HER4","HSP-27","IGF1-R","ILK","INH","MDM2","MIB-1","O13","P16","P21","P27","P53","P57","P63","P-AKT","PGI","PR","PSA","AR","CK19","COX2","Galectin-3","HBME-1","MLH1","MSH2","MSH6","P504S","P75-NTR","PMS2","TOPO-II","TS106","RET","S100","Syntrophin","TDT","TSH","TTF-1","TG","UPA","VEGF","WT1")
all_col_names=c("AAT","AMF-R","Aurora-A","Aurora-B","Aurora-C","Bcl-2","CTNNB1","CAIX","CAV-1","CDX2","c-kit","CR3","Cyclin-D1","Cyclin-E","Caveolin","Clusterin","E-CAD","EGFR","ER","HER2","HER3","HER4","HSP-27","IGF1-R","ILK","INH","MDM2","MIB-1","O13","P16","P21","P27","P53","P57","P63","P-AKT","PGI","PR","PSA","AR","CK19","COX2","Galectin-3","HBME-1","MLH1","MSH2","MSH6","P504S","P75-NTR","PMS2","TOPO-II","TS106","RET","S100","Syntrophin","TDT","TSH","TTF-1","TG","UPA-R","VEGF","WT1")

x=as.matrix(marker_data)
#rownames(x)=pathology_patient
rownames(x)=pathology_patient_new
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=all_col_names, col=score_colors, RowSideColors=patientcolors, cexRow=0.8, cexCol=0.60)
dev.off()

#Create a vector assigning color to each marker based on whether it is up- or down-regulated for a second color side bar
sigmarkerchange=c("DOWN","DOWN","UP","DOWN","UP","UP","DOWN","DOWN")
sigmarkercolors=sigmarkerchange
sigmarkercolors[sigmarkercolors=="DOWN"]="#FFA500" #Orange
sigmarkercolors[sigmarkercolors=="UP"]="#FFFF00" #Yellow

#Now get matrix for just markers that were sig associated with cancer status
#pdf("ATC_all_deconvoluted_data_62markers_23JAN07_matched_heatmap_wkey_MHsig.pdf")
pdf("ATC_all_deconvoluted_data_62markers_23JAN07_matched_heatmap_wkey_MHsig.2.pdf")
#score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("BCL2","B.Catenin","MIB1","ECAD","P53","TOPO2","VEGF","Thyroglobulin")]
#Specify column names separately
#sig_col_names=c("BCL2","B-Catenin","MIB1","ECAD","P53","TOPO2","VEGF","Thyroglobulin")
sig_col_names=c("Bcl-2","CTNNB1","MIB-1","E-CAD","P53","TOPO-II","VEGF","TG")
x2=as.matrix(sig_marker_data)
#rownames(x2)=pathology_patient
rownames(x2)=pathology_patient_new
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, col=score_colors, ColSideColors=sigmarkercolors, RowSideColors=patientcolors, cexRow=1.2, cexCol=1.2)
dev.off()


#######################################
#Unused code
#for problems with labels not fitting on page, try using par() with ps= and cex= options.

#For creating a normal hierarchical tree
#distances=dist(marker_data, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#hc=hclust(distances, method = "complete")
#plot(hc,labels=pathology)

