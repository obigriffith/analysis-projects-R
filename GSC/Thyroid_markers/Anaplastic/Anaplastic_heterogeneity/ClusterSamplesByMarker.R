library("gplots") #For advanced heatmaps
library(RColorBrewer)

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/Anaplastic_heterogeneity")

data=read.table("ATC_all_deconvoluted_data_52markers_25OCT06_spindle_epith_giant_matched.txt", header = TRUE, na.strings = "NA", sep="\t")

marker_data=data[,7:57]


pathology=as.vector(data[,6])
pathology_patient=paste(pathology,data[,1], sep="_")

#Create a vector assigning a color to each patient based on pathology for the color side bar
#Use brewer.pal to obtain color palettes. Use display.brewer.all() for selection
#colors=brewer.pal(9,"Set1")
colors=brewer.pal(8,"Dark2")
patientcolors=pathology
patientcolors[patientcolors=="spindled"]=colors[1]
patientcolors[patientcolors=="giant_cell"]=colors[2]
patientcolors[patientcolors=="epithelioid"]=colors[3]

#Alternatively, assign a color for primary vs secondary foci
foci=data[,2]
focicolors=foci
focicolors[focicolors==1]=colors[1]
focicolors[focicolors==2]=colors[2]

#Assign colors for the marker scores
score_col_pal=brewer.pal(9, "Blues")
#Choose colors for the five possible scores (0,1,2,3,4)
#score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#score_colors=c(score_col_pal[2],score_col_pal[4],score_col_pal[6],score_col_pal[8],score_col_pal[9]) 
score_colors=c(score_col_pal[1],score_col_pal[3],score_col_pal[5],score_col_pal[7],score_col_pal[9]) 

#Set marker names to be consistent with previous publications
#col_names=colnames(marker_data)
col_names=c("AAT","AMF-R","Aurora-A","Aurora-B","Aurora-C","Bcl-2","CTNNB1","CAIX","CDX2","c-kit","CR3","Cyclin-D1","Cyclin-E","Clusterin","E-CAD","ER","EGFR","HER2","HER3","HER4","HSP-27","IGF1-R","ILK","INH","MDM2","MIB-1","O13","P16","P21","P27","P53","P57","P63","PR","PSA","AR","CK19","Galectin-3","P75-NTR","PMS2","TOPO-II","TS106","RET","S100","TDT","TSH","TTF-1","UPA-R","VEGF","WT1","TG")

#Simplify patient names for publication.
#row_names=pathology_patient
row_names=c("P1","S1","P2","S2","P3","S3","P4","S4","P5","S5","P6","S6")
row_names2=c("spindled_1","giant_cell_1","epithelioid_2","spindled_2","spindled_3","giant_cell_3","spindled_4","epithelioid_4","spindled_5","epithelioid_5","spindled_6","epithelioid_6")

#To create a heatmap of all data
#Try different color scheme and add a score key
pdf("ATC_52markers_25OCT06_heterogeneity_heatmap_by_type_wkey.pdf")
x=as.matrix(marker_data)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", col=score_colors, RowSideColors=patientcolors, labCol=col_names, labRow=row_names, cexRow=1, cexCol=0.8)
dev.off()

#Try with different sidebar 
pdf("ATC_52markers_25OCT06_heterogeneity_heatmap_by_Foci_wkey.pdf")
x=as.matrix(marker_data)
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", col=score_colors, RowSideColors=focicolors, labCol=col_names, labRow=row_names2, cexRow=0.8, cexCol=0.8)
dev.off()
