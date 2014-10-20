library("gplots")
library("heatmap.plus")

#Gene
#datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/33libs/matrix/Matrix_GeneExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/33libs/matrix/Expressed_GeneExpression_v53.txt"

#Transcript
datafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/33libs/matrix/Matrix_TranscriptExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/33libs/matrix/Expressed_TranscriptExpression_v53.txt"

outdir="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/33libs/"
#outfile = "BCCL_Gene_expression_33libs_heatmap.pdf"
outfile = "BCCL_Transcript_expression_33libs_heatmap.pdf"
#outfile = "BCCL_Gene_expression_33libs_heatmap_HER2.pdf"
#outfile = "BCCL_Transcript_expression_33libs_heatmap_HER2.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
exp_status=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
setwd(outdir)
header=colnames(raw_data)

pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 2 #Minimum coefficient of variation (0.7 recommended?)
cov_max = 10 #Maximum cov

#All libraries
#libs=header[5:37]
#lib_names=c("184A1","184B5","600MPE","BT474","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC3153","HCC70","M4A4","MCF10A","MCF10F","MCF7","MDAMB134VI","MDAMB157","MDAMB231","Normal","SKBR3","SUM149PT","SUM52PE","T47D","ZR7530","ZR75B")
#lib_types=c("Basal_NM","Basal_NM","Luminal","Luminal","ClaudinLow","Luminal","Basal","Failed","ClaudinLow","Luminal","Luminal","Basal","Basal","Basal","Basal","Basal","Luminal","Basal","Basal","Other","Basal_NM","Basal_NM","Luminal","Luminal","ClaudinLow","ClaudinLow","Normal","Luminal","Basal","Luminal","Luminal","Luminal","Luminal")
#her2_status=c("HER2NEG","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2POS","HER2POS","HER2NEG","HER2NEG","Unknown","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","Normal","HER2POS","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG")

#All libraries except normal, failed, other
libs=header[c(5:11,13:23,25:30,32:37)]
lib_names=c("184A1","184B5","600MPE","BT474","BT549","CAMA1","HCC1143","HCC1395","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC3153","HCC70","MCF10A","MCF10F","MCF7","MDAMB134VI","MDAMB157","MDAMB231","SKBR3","SUM149PT","SUM52PE","T47D","ZR7530","ZR75B")
lib_types=c("Basal_NM","Basal_NM","Luminal","Luminal","ClaudinLow","Luminal","Basal","ClaudinLow","Luminal","Luminal","Basal","Basal","Basal","Basal","Basal","Luminal","Basal","Basal","Basal_NM","Basal_NM","Luminal","Luminal","ClaudinLow","ClaudinLow","Luminal","Basal","Luminal","Luminal","Luminal","Luminal")
her2_status=c("HER2NEG","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2POS","HER2POS","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2POS","HER2NEG","HER2NEG","HER2NEG","HER2NEG","HER2NEG")

#Set up colors for sidebars
colors=rainbow(8)
#Subtype
subtype_colors=lib_types
subtype_colors[subtype_colors=="Basal"]=colors[1]
subtype_colors[subtype_colors=="Luminal"]=colors[3]
subtype_colors[subtype_colors=="ClaudinLow"]=colors[5]
subtype_colors[subtype_colors=="Basal_NM"]=colors[7]
subtype_colors[subtype_colors=="Normal"]=colors[6]
subtype_colors[subtype_colors=="Other"]=colors[4]
subtype_colors[subtype_colors=="Failed"]=colors[4]

#HER2 status
her2_colors=her2_status
her2_colors[her2_colors=="HER2NEG"]=colors[2]
her2_colors[her2_colors=="HER2POS"]=colors[8]
her2_colors[her2_colors=="Normal"]=colors[6]
her2_colors[her2_colors=="Unknown"]=colors[4]

#Define a percent expressed function and filter out features with less than minimum
w=exp_status[,libs]
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}
pe_data=apply(w, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data=raw_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
y=data[,libs]
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}
cov_data=apply(y, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data=data[passed_cov,]

#Grab just library data for genes remaining after filtering
z=data[,libs]

#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(z+1)

#Heatmap
pdf(file=outfile)
x=t(as.matrix(z))
#x=t(x)
#main_title="RNAseq gene-level expression"
main_title="RNAseq transcript-level expression"
rlab=cbind(subtype_colors,her2_colors)
colnames(rlab) = c("Subtype","Her2")
heatmap.plus(x, na.rm = TRUE, scale="none", main=main_title, RowSideColors=rlab, labRow=lib_names, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
dev.off()


#heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, RowSideColors=subtype_colors, labRow=lib_names, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
#heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, RowSideColors=her2_colors, labRow=lib_names, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
#Other color possibilities: col=rich.colors(75, palette="blues")
