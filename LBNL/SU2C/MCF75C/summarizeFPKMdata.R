library("gplots")

#Transcript-level fpkm data from cufflinks
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/MCF75C/TophatCufflinks/formattedFiles/MCF75C_fpkm.txt"

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/MCF75C/TophatCufflinks/plotsTables"


heatmap_outfile = "FPKM_expression_all_libs_heatmap.pdf"

setwd(outdir)


#Parameters
pe_thresh = 0.15 #Minimum percent libraries "expressed"
cov_min = 2 #Minimum coefficient of variation (0.7 recommended?)
#cov_min = 0 #Minimum coefficient of variation (0.7 recommended?)
cov_max = 10 #Maximum cov

#Import data
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", as.is=c(1:2))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:2]
raw_data=raw_data_import[,3:length(colnames(raw_data_import))]
libnames=colnames(raw_data)
libnames_clean=c("Control","E2","4OHT","E2+4OHT","PP2","E2+PP2")

#Create a corresponding dataframe with expression status
exp_data=raw_data
exp_data[exp_data>=1]=1
exp_data[exp_data<1]=0


#Set up colors for sidebars, etc
colors=rainbow(4)
#Model
model_colors=c(rep(colors[1],8),rep(colors[2],7),rep(colors[3],7),rep(colors[4],7))


#Define a percent expressed function and filter out features with less than minimum
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}

pe_data=apply(exp_data, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data_filt=raw_data[passed_pe,]
feat_data_filt=raw_feat_data[passed_pe,]

cov_data=apply(data_filt, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data_filt=data_filt[passed_cov,]
feat_data_filt=feat_data_filt[passed_cov,]

gene_names=feat_data_filt[,2]

#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(data_filt+1)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}


x=t(as.matrix(z))
main_title="Cufflinks transcript-level expression"
par(cex.main=.7)
pdf(file=heatmap_outfile)
#heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=gene_names, cexRow=0.8, cexCol=0.5, col=rev(heat.colors(75)))
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, labRow=libnames_clean, labCol=FALSE, cexRow=0.8, cexCol=0.5, col=rev(heat.colors(75)))
dev.off()

