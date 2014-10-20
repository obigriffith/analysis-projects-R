library("gplots") #For advanced heatmaps

setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/Anaplastic/complete_cohort_analysis")
datafile="ATC_all_deconvoluted_data_62markers_23JAN07_nocalcpos_primfoci.txt"

data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

type=data[,7]
marker_data=data[,8:69]

diff_foci=data[,6]
diff_foci[diff_foci==0]="N"
diff_foci[diff_foci==1]="Y"

#Create a vector assigning a color to each patient based on pathology for the color side bar
typecolors=as.vector(type)
typecolors[typecolors=="epithelioid"]="#9ACD32" #yellowgreen
typecolors[typecolors=="spindled"]="#006400" #darkgreen
typecolors[typecolors=="squamoid"]="#FF0000" #red

diff_foci_colors=diff_foci
diff_foci_colors[diff_foci_colors=="N"]="#9ACD32" #yellowgreen
diff_foci_colors[diff_foci_colors=="Y"]="#006400" #darkgreen


###To create a heatmap of all data###
#Try different color scheme and add a score key
#Choose colors for the five possible scores (0,1,2,3,4)
score_colors=c("#F0F8FF","#B9D3EE","#00BFFF","#483D8B","#00008B") 
#Specify column names
all_col_names=colnames(marker_data)

pdf("ATC_all_deconvoluted_data_62markers_23JAN07_nocalcpos_primfoci_heatmap_wkey.pdf")
x=as.matrix(marker_data)
rownames(x)=type
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=typecolors, cexRow=0.8, cexCol=0.60)
dev.off()

pdf("ATC_all_deconvoluted_data_62markers_23JAN07_nocalcpos_primfoci_heatmap_wkey_df.pdf")
x=as.matrix(marker_data)
rownames(x)=type
heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=FALSE, labCol=all_col_names, col=score_colors, RowSideColors=diff_foci_colors, cexRow=0.8, cexCol=0.60)
dev.off()
