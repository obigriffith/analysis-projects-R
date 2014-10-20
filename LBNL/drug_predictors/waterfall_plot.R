#Set file names and working dir
#drug="BIBW2992"
drug="GSK_Tykerb" #Lapatinib

outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/temp/"
dir.create(outdir) #Creates directory if not already in existence
setwd(outdir)
waterfallfile="waterfall.pdf"

#Import cell line data
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Filter down to core set of cell lines - must have drug data and at least one other molecular profiling data type
core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")
cell_line_data=raw_cell_line_import[core_cell_lines,]

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

#For drug of interest divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)

#Get sorted data and list of cell lines for waterfall plot
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]

#Determine mean GI50 to show threshold of waterfall plot
mean_cutoff=mean(drug_data_interest, na.rm=TRUE)
resistants=which(drug_data_interest_sorted<=mean_cutoff)
sensitives=which(drug_data_interest_sorted>mean_cutoff)
num_sensitive=length(sensitives)

#Set colors for subtype color sidebar
subtypes=as.vector(cell_line_data[cell_lines_sorted,"BCCLclassification2"])
subtype_colors=subtypes
subtype_colors[subtype_colors=="ERBB2Amp"]="blue"
subtype_colors[subtype_colors=="Basal"]="red"
subtype_colors[subtype_colors=="Claudin-low"]="green"
subtype_colors[subtype_colors=="Luminal"]="black"

#Create corresponding histogram/waterfall plot for drug response data
ymin=floor(min(drug_data_interest_sorted))
ymax=ceiling(max(drug_data_interest_sorted))
barplot_values=barplot(drug_data_interest_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
pdf(file=waterfallfile)
barplot(drug_data_interest_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, col=subtype_colors, ylab="-log10(GI50)", main=drug)
#draw line midway between bars straddling sensitive/resistant threshold
abline(v=barplot_values[num_sensitive]+(barplot_values[num_sensitive+1]-barplot_values[num_sensitive])/2, lwd=2, col="black")
legend("topright", legend=c("Luminal","Basal","Claudin-low","ERBB2Amp"), fill=c("black","red","green","blue"), cex=0.9, bty="n")
box()
dev.off()
