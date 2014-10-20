library("gplots") #For advanced heatmaps

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis/data_files")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"

data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Keep only rows where pathology is 0, 1, 2, or 3. Note we are only excluding M cases not HCC as previously done
pathology=data[,13]
data=data[pathology<4,]
pathology=data[,13] #reset pathology after exclusion
marker_data=data[,36:92] #ungrouped values

pathology[pathology==0]="B"
pathology[pathology==1]="F"
pathology[pathology==2]="P"
pathology[pathology==3]="H"
pathology[pathology==4]="M"

#Create a vector assigning a color to each patient based on pathology for the color side bar
patientcolors=pathology
patientcolors[patientcolors=="B"]="#9ACD32" #yellowgreen
#patientcolors[patientcolors=="F"]="#FF0000" #red
patientcolors[patientcolors=="F"]="#006400" #darkgreen
patientcolors[patientcolors=="P"]="#006400" #darkgreen
patientcolors[patientcolors=="H"]="#006400" #darkgreen

###To create a heatmap of all data###
#Try different color scheme and add a score key
pdf("BenignMalignant_58markers_23JAN08_heatmap_wkey.pdf")
#Choose colors for the five possible scores (0,1,2,3,4)
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names separately
all_col_names=c("AAT","AMF-R","AR","Aurora-A","Aurora-C","Bcl-2","CTNNB1","CDX2","CK19","c-kit","COX2","CR3","Cyclin-D1","Cyclin-E","Caveolin","CAV-1","Clusterin","E-CAD","ER","Galectin-3","HBME-1","EGFR","HER2","HER3","HER4","HSP-27","IGFBP2","IGFBP5","INH","MIB-1","MDM2","MLH1","MRAS","O13","P16","P21","P27","P53","P57","P63","P75-NTR","PGI","PMS2","PR","PSA","RET","S100","Syntrophin","TDT","TG","TOPO-II","TS106","TSH","TTF-1","VEGF","WT1","P504S")
#all_col_names=c("AAT","AMF-R","AR","Aurora-A","Aurora-C","Bcl-2","CTNNB1","CDX2","CK19","c-kit","COX2","CR3","Cyclin-D1","Cyclin-E","Caveolin","CAV-1","Clusterin","E-CAD","ER","Galectin-3","HBME-1","EGFR","HER2","HER3","HER4","HSP-27","IGFBP2","IGFBP5","INH","KI67","MDM2","MLH1","MRAS","O13","P16","P21","P27","P53","P57","P63","P75-NTR","PGI","PMS2","PR","PSA","RET","S100","Syntrophin","TDT","TG","TOPO-II","TS106","TSH","TTF-1","VEGF","WT1","P504S")
#all_col_names=colnames(marker_data)

x=as.matrix(marker_data)
rownames(x)=pathology
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=patientcolors, cexRow=0.8, cexCol=0.60)
dev.off()

###All significant markers###
#Now get matrix for just markers that were sig associated with cancer status (significant in any test: Chi-square group1, group2 or ungrouped Mann-Whitney)
#Create a vector assigning color to each marker based on whether it is up- or down-regulated for a second color side bar
sigmarkerchange=c("UP","UP","DOWN","UP","UP","UP","UP","DOWN","UP","UP","UP","DOWN","UP","UP","UP","UP","UP","UP","UP","DOWN","UP","UP","UP","UP","DOWN","DOWN","UP","UP","UP","DOWN","UP","UP","UP","DOWN","UP")
sigmarkercolors=sigmarkerchange
sigmarkercolors[sigmarkercolors=="DOWN"]="#FFA500" #Orange
sigmarkercolors[sigmarkercolors=="UP"]="#FFFF00" #Yellow

pdf("BenignMalignant_58markers_23JAN08_heatmap_wkey_sig.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME","BCL2","CYCLIND1","Caveolin1","CYCLINE","ECAD","CR3","Clusterin","IGFBP5","P21","IGFBP2","BetaCatenin","HER4","TG","KI67","Caveolin","AuroraC","S100","MRAS","CKIT","HER3","RET","AMFR","MLH1","AAT","TTF1","PGI","HSP27","Syntrophin")]
#Specify column names separately
sig_col_names=c("Galectin-3","CK19","VEGF","Aurora-A","P16","AR","HBME-1","Bcl-2","Cyclin-D1","CAV-1","Cyclin-E","E-CAD","CR3","Clusterin","IGFBP5","P21","IGFBP2","CTNNB1","HER4","TG","MIB-1","Caveolin","Aurora-C","S100","MRAS","c-kit","HER3","RET","AMF-R","MLH1","AAT","TTF-1","PGI","HSP-27","Syntrophin")
x2=as.matrix(sig_marker_data)
rownames(x2)=pathology
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, ColSideColors=sigmarkercolors, RowSideColors=patientcolors, cexRow=0.9, cexCol=0.9)
dev.off()

###Top 10 Most significant markers from Mann-Whitney###
#Create a vector assigning color to each marker based on whether it is up- or down-regulated for a second color side bar
sigmarkerchange=c("UP","UP","DOWN","UP","UP","UP","UP","DOWN","UP","UP")
sigmarkercolors=sigmarkerchange
sigmarkercolors[sigmarkercolors=="DOWN"]="#FFA500" #Orange
sigmarkercolors[sigmarkercolors=="UP"]="#FFFF00" #Yellow

pdf("BenignMalignant_58markers_23JAN08_heatmap_wkey_top10_mostsig.pdf")
#score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME","BCL2","CYCLIND1","Caveolin1")]
#Specify column names separately
sig_col_names=c("Galectin-3","CK19","VEGF","Aurora-A","P16","AR","HBME-1","Bcl-2","Cyclin-D1","CAV-1")
x2=as.matrix(sig_marker_data)
rownames(x2)=pathology
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, ColSideColors=sigmarkercolors, RowSideColors=patientcolors, cexRow=1, cexCol=1)
dev.off()


###Top 7 Most significant markers from Mann-Whitney###
#Now get matrix for just markers that were really significant
#Create a vector assigning color to each marker based on whether it is up- or down-regulated for a second color side bar
sigmarkerchange=c("UP","UP","DOWN","UP","UP","UP","UP")
sigmarkercolors=sigmarkerchange
sigmarkercolors[sigmarkercolors=="DOWN"]="#FFA500" #Orange
sigmarkercolors[sigmarkercolors=="UP"]="#FFFF00" #Yellow

pdf("BenignMalignant_58markers_23JAN08_heatmap_wkey_mostsig.pdf")
#score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME")]
#Specify column names separately
sig_col_names=c("Galectin-3","CK19","VEGF","Aurora-A","P16","AR","HBME-1")
x2=as.matrix(sig_marker_data)
rownames(x2)=pathology
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, ColSideColors=sigmarkercolors, RowSideColors=patientcolors, cexRow=1.2, cexCol=1.05)
dev.off()

#######################################
#Unused code
#for problems with labels not fitting on page, try using par() with ps= and cex= options.

#For creating a normal hierarchical tree
#distances=dist(marker_data, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#hc=hclust(distances, method = "complete")
#plot(hc,labels=pathology)

