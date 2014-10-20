library("gplots") #For advanced heatmaps

#Set working directory and filenames for Input/output
#setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimfoci_vs_DTC")
#setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCdiff_Foci_vs_DTC")
#setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimFociNoDiff_vs_DTC")
setwd("C:/Users/Obi/Documents/GSC/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis/ATC_vs_DTC/ATCprimFociNoDiff_vs_DTC")

#datafile="ATCprimfoci_vs_DTC.txt"
datafile="ATCprimFociNoDiff_vs_DTC.txt"
#datafile="ATCdiff_Foci_vs_DTC.txt"

data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Get diagnostic target variable and specify as factor/categorical
type=data[,7]
type[type==0]="ATC"
type[type==1]="DTC"
target=as.factor(type)

#Get marker data
marker_data=data[,10:63] #ungrouped values

#Create a vector assigning a color to each patient based on pathology for the color side bar
typecolors=as.vector(type)
typecolors[typecolors=="ATC"]="#9ACD32" #yellowgreen
typecolors[typecolors=="DTC"]="#006400" #darkgreen

all_marker_data=marker_data[,c("AAT","AMF","AR","AuroraA","AuroraC","BCL2","BetaCatenin","CDX2","CK19","CKIT","COX2","CR3","CYCLIND1","CYCLINE","Caveolin","CAV1","Clusterin","ECAD","ER","Galectin","HBME1","EGFR","HER2","HER3","HER4","HSP27","INH","MIB1","MDM2","MLH1","O13","P16","P21","P27","P53","P57","P63","P75_NTR","PGI","PMS2","PR","PSA","RET","S100","Syntrophin","TDT","Thyroglobulin","TOPO2","TS106","TSH","TTF1","VEGF","WT1","P504S")]

###Heatmaps### MAKE SURE YOU HAVE THE CORRECT DATA LOADED FROM ABOVE BEFORE RUNNING EACH HEATMAP BELOW
#ATCprimfoci_vs_DTC, All markers
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names
all_col_names=colnames(marker_data)
pdf("ATCprimfoci_vs_DTC_heatmap_wkey.pdf")
x=as.matrix(all_marker_data)
rownames(x)=type
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=typecolors, cexRow=0.8, cexCol=0.60)
dev.off()

#ATCprimfoci_vs_DTC, sig markers
pdf("ATCprimfoci_vs_DTC_heatmap_wkey_sig.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
sig_marker_data=marker_data[,c("AAT","AMF","AR","BetaCatenin","CK19","COX2","CR3","CYCLINE","Clusterin","ECAD","HBME1","HER3","HER4","MIB1","P21","P27","P53","P63","PGI","PR","RET","Syntrophin","Thyroglobulin","TOPO2","TS106","TTF1")]
#Specify column names separately
sig_col_names=c("AAT","AMF-R","AR","CTNNB1","CK19","COX2","CR3","Cyclin-E","Clusterin","E-CAD","HBME-1","HER3","HER4","MIB-1","P21","P27","P53","P63","PGI","PR","RET","Syntrophin","TG","TOPO-II","TS106","TTF-1")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()

#ATCprimfoci_vs_DTC, Wiseman top7 markers (Benign/Malignant diagnostic panel)
pdf("ATCprimfoci_vs_DTC_heatmap_wkey_Wisemantop7.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")]
#Specify column names separately
sig_col_names=c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()


#################
#ATCdiff_Foci_vs_DTC, All markers
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names
all_col_names=colnames(marker_data)
pdf("ATCdiff_Foci_vs_DTC_heatmap_wkey.pdf")
x=as.matrix(all_marker_data)
rownames(x)=type
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=typecolors, cexRow=0.8, cexCol=0.60)
dev.off()

#ATCdiff_Foci_vs_DTC, sig markers
pdf("ATCdiff_Foci_vs_DTC_heatmap_wkey_sig.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("BetaCatenin","TOPO2","TTF1","CKIT","CK19","VEGF","BCL2","CR3","Galectin","P27","HBME1","PGI","CAV1","TS106","WT1","AuroraA","Clusterin","RET","MDM2","EGFR","Syntrophin")]
#Specify column names separately
sig_col_names=c("BetaCatenin","TOPO2","TTF1","CKIT","CK19","VEGF","BCL2","CR3","Galectin","P27","HBME1","PGI","CAV1","TS106","WT1","AuroraA","Clusterin","RET","MDM2","EGFR","Syntrophin")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()

#ATCdiff_Foci_vs_DTC, Wiseman top7 markers (Benign/Malignant diagnostic panel)
pdf("ATCdiff_Foci_vs_DTC_heatmap_wkey_Wisemantop7.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")]
#Specify column names separately
sig_col_names=c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()



#################
#ATCprimFociNoDiff_vs_DTC, All markers
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names
all_col_names=colnames(marker_data)
pdf("ATCprimFociNoDiff_vs_DTC_heatmap_wkey.pdf")
x=as.matrix(all_marker_data)
rownames(x)=type
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=typecolors, cexRow=0.8, cexCol=0.60)
dev.off()

#ATCprimFociNoDiff_vs_DTC, sig markers
pdf("ATCprimFociNoDiff_vs_DTC_heatmap_wkey_sig.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
sig_marker_data=marker_data[,c("AAT","AMF","AR","BetaCatenin","Caveolin","CK19","COX2","CR3","CYCLINE","Clusterin","ECAD","HBME1","EGFR","HER3","HER4","MIB1","P21","P53","RET","S100","Syntrophin","Thyroglobulin","TOPO2","TS106","TTF1")]
#Specify column names separately
sig_col_names=c("AAT","AMF-R","AR","CTNNB1","Caveolin","CK19","COX2","CR3","Cyclin-E","Clusterin","E-CAD","HBME-1","EGFR","HER3","HER4","MIB-1","P21","P53","RET","S100","Syntrophin","TG","TOPO-II","TS106","TTF-1")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()

#ATCprimFociNoDiff_vs_DTC, Wiseman top7 markers (Benign/Malignant diagnostic panel)
pdf("ATCprimFociNoDiff_vs_DTC_heatmap_wkey_Wisemantop7.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")]
#Specify column names separately
sig_col_names=c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME1")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()

#ATCprimFociNoDiff_vs_DTC, Wiseman top8 markers (ATC transformation paper)
pdf("ATCprimFociNoDiff_vs_DTC_heatmap_wkey_Wiseman_trans_top8.pdf")
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B") 
sig_marker_data=marker_data[,c("BCL2","BetaCatenin","ECAD","VEGF","MIB1","P53","Thyroglobulin","TOPO2")]
#Specify column names separately
sig_col_names=c("Bcl-2","CTNNB1","E-CAD","VEGF","MIB-1","P53","TG","TOPO-II")
x2=as.matrix(sig_marker_data)
rownames(x2)=type
heatmap.2(x2, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=sig_col_names, labRow=FALSE, col=score_colors, RowSideColors=typecolors, cexRow=0.9, cexCol=0.9)
dev.off()

