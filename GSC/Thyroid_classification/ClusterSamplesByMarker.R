library("gplots") #For advanced heatmaps

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_classification/clustering/BenignMalignant")

data=read.table("BenignMalignant_35markers_2AUG06_noHCCorM.txt", header = TRUE, na.strings = "NA", sep="\t")

#Skip col 44 (Galectin) and col 69 (EGFRBlock)
marker_data=data[,c(35:43,45:68)] #Skip col 44 (Galectin) and col 69 (EGFRBlock)

pathology=data[,12]
pathology[pathology==0]="B"
pathology[pathology==1]="F"
pathology[pathology==2]="P"

#Create a vector assigning a color to each patient based on pathology for the color side bar
patientcolors=pathology
patientcolors[patientcolors=="B"]="#FFFF00" #yellow
patientcolors[patientcolors=="F"]="#FF0000" #red
patientcolors[patientcolors=="P"]="#00FF00" #green


#To create a heatmap of all data
pdf("BenignMalignant_33markers_2AUG06_noHCCorM_heatmap_wkey.pdf")
#Choose colors for the five possible scores (0,1,2,3,4)
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
x=as.matrix(marker_data)
rownames(x)=pathology
heatmap.2(t(x), na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", col=score_colors, ColSideColors=patientcolors, labCol=FALSE, cexRow=0.8)
dev.off()

#Now get matrix for just markers that were significantly associated with cancer status
pdf("BenignMalignant_33markers_2AUG06_noHCCorM_heatmap_wkey_sigmarkers_only.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
sig_marker_data=marker_data[,c("AR","VEGF","P16","BCL2block","CYCLIND1","CYCLINE","P21","S100","KI67","CKIT","HER3","AMFR","AAT","HER1")]
x2=as.matrix(sig_marker_data)
rownames(x2)=pathology
heatmap.2(t(x2), na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", col=score_colors, ColSideColors=patientcolors, labCol=FALSE, cexRow=0.8)
dev.off()

#Now get matrix for just the top 5 best markers that were significantly associated with cancer status
pdf("BenignMalignant_33markers_2AUG06_noHCCorM_heatmap_wkey_bestsigmarkers_only.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") #Without AAT there are no markers that have score=4
sig_marker_data=marker_data[,c("AR","VEGF","P16","BCL2block","CYCLIND1")]
x3=as.matrix(sig_marker_data)
rownames(x3)=pathology
heatmap.2(t(x3), na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", col=score_colors, ColSideColors=patientcolors, labCol=FALSE, cexRow=0.8)
dev.off()

#######################################
#Unused code
#for problems with labels not fitting on page, try using par() with ps= and cex= options.

#For creating a normal hierarchical tree
#distances=dist(marker_data, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#hc=hclust(distances, method = "complete")
#plot(hc,labels=pathology)

