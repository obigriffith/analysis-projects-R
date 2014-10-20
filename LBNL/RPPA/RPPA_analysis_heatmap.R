library("gplots") #For advanced heatmaps

setwd("C:/Users/Obi/Documents/Projects/HER2_Breast_Cancer/RPPA/R_analysis")

#Load data
datafile="JoeGray480_MedianCentered_179Ab.clean.data.txt"
rawdata=read.table(datafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t")

#Load protein details
protein_metafile="JoeGray480_MedianCentered_179Ab.clean.protein_meta.txt"
protein_meta=read.table(protein_metafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t", as.is=TRUE)

#Load cell line details
lines_metafile="JoeGray480_MedianCentered_179Ab.clean.lines_meta.txt"
lines_meta_raw=read.table(lines_metafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t")

##choose either (1) or (2)
#1. All lines
#data=rawdata
#lines_meta=lines_meta_raw

#2. HER2+ lines only (Amplified only)
data=rawdata[rawdata[,"HER2_status"]=="Amplified",]
lines_meta=lines_meta_raw[lines_meta_raw[,"HER2_status"]=="Amplified",]

#choose either (A) or (B)
#Define comparison groups
#A. Synergy vs non-synergy
groupA_lines=which(lines_meta[,"Lap_AKT_Synergy"]=="Y")
groupB_lines=which(lines_meta[,"Lap_AKT_Synergy"]=="N")

#B. mutant (PI3K or PTEN) vs wildtype
#groupA_lines=which(lines_meta[,"PI3K_OR_PTEN"]=="mutant")
#groupB_lines=which(lines_meta[,"PI3K_OR_PTEN"]=="wt")

##Specify output files
#Choose only one of the following file names (1A, 1B, 2A, or 2B, etc)
#1A. All lines, synergy vs non-synergy 
#pdf(file="synergy_vs_nonsynergy_all.heatmap.pdf")

#1B. All lines, mutant vs wild-type 
#pdf(file="mutant_vs_wildtype_all.heatmap.pdf")

#2A. HER2+ lines only (Amplified only), synergy vs non-synergy
file1="synergy_vs_nonsynergy_HER2amp.heatmap.ALL.pdf"
file2="synergy_vs_nonsynergy_HER2amp.heatmap.highCOV.pdf"
file3="synergy_vs_nonsynergy_HER2amp.heatmap.lowCOV.pdf"
file4="synergy_vs_nonsynergy_HER2amp.dists.pdf"

#2B. HER2+ lines only (Amplified only), mutant vs wild-type
#pdf(file="mutant_vs_wildtype_HER2amp.heatmap.pdf")

#Specify some variables for filtering purposes
treatments=c("DMSO","Lapatinib","AKTi","Lap_AKTi")
timepoints=c(0.5,1.0,2.0,4.0,8.0,24.0,48.0,72.0)
lines=rownames(lines_meta)

protein_names=protein_meta[1,]
RPPA_qualities=protein_meta[3,]
RPPA_qual_col=protein_meta[2,]
RPPA_qual_col[RPPA_qual_col=="F"]="red"
RPPA_qual_col[RPPA_qual_col=="M"]="red"
RPPA_qual_col[RPPA_qual_col=="P"]="black"

#Separate abundance values from first 4 descriptor columns
descriptor=data[,1:9]
data=data[,10:188]

proteins=colnames(data)

#Now, average across cell lines in each comparison group and examine across each time point individually
#Compile list of treatment/timepoint combinations
treatment_time=vector(length=length(treatments)*length(timepoints))
n=0
for (i in 1:length(treatments)){
 for (j in 1:length(timepoints)){
  n=n+1
  treatment_time[n]=paste(treatments[i],timepoints[j],sep="_")
 }
}
treatment_time_array = array(0, dimnames = list(treatment_time, c("treatment","time")), dim=c(length(treatment_time),2))


#Mean protein abundance across all cell lines (for each category) for each treatment/time
groupA_mean_data = array(0, dimnames = list(treatment_time, colnames(data)), dim=c(length(treatment_time),length(colnames(data))))
groupB_mean_data = array(0, dimnames = list(treatment_time, colnames(data)), dim=c(length(treatment_time),length(colnames(data))))

n=0
for (i in 1:length(treatments)){
 for (j in 1:length(timepoints)){
  n=n+1
  treatment_time_data=data[descriptor[,"time_hrs"]==timepoints[j] & descriptor[,"treatment"]==treatments[i],]
  #Stratify by group
  groupA_treatment_time_data=treatment_time_data[groupA_lines,]
  groupB_treatment_time_data=treatment_time_data[groupB_lines,]
  mean_groupA_treatment_time_data=mean(groupA_treatment_time_data)
  mean_groupB_treatment_time_data=mean(groupB_treatment_time_data)
  groupA_mean_data[n,]=mean_groupA_treatment_time_data
  groupB_mean_data[n,]=mean_groupB_treatment_time_data
  treatment_time_array[n,]=c(treatments[i],timepoints[j])
 }
}

#Calculate difference between means for groupA (e.g., synergy) and groupB (e.g., non-synergy)
groupAB_mean_diff=groupA_mean_data-groupB_mean_data

#Cluster and create heatmap
colors=rainbow(4)
#Treatment
treatment_colors=treatment_time_array[,"treatment"]
treatment_colors[treatment_colors=="DMSO"]=colors[1]
treatment_colors[treatment_colors=="Lapatinib"]=colors[2]
treatment_colors[treatment_colors=="AKTi"]=colors[3]
treatment_colors[treatment_colors=="Lap_AKTi"]=colors[4]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_min = 0.1 #Chosen by experimentation
y=groupAB_mean_diff+abs(min(groupAB_mean_diff)) #Scale up so that smallest value is zero
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}
cov_data=apply(y, 2, cov_fun)



#All proteins
pdf(file=file1)
main_title="Mean diff (syn - nonsyn): All"
x=t(as.matrix(groupAB_mean_diff))
heatmap_image=heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="row", main=main_title, ColSideColors=treatment_colors, cexRow=0.3, cexCol=0.8, col=bluered, symbreaks=TRUE)
dev.off()

#High COV proteins
passed_cov = which(cov_data >= cov_min)
groupAB_mean_diff_filt=groupAB_mean_diff[,passed_cov]
pdf(file=file2)
main_title="Mean diff (syn - nonsyn): COV > 0.1"
x=t(as.matrix(groupAB_mean_diff_filt))
heatmap_image=heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="row", main=main_title, ColSideColors=treatment_colors, cexRow=0.7, cexCol=0.8, col=bluered, symbreaks=TRUE)
dev.off()

#Low COV proteins
passed_cov = which(cov_data < cov_min)
groupAB_mean_diff_filt=groupAB_mean_diff[,passed_cov]
pdf(file=file3)
main_title="Mean diff (syn - nonsyn): COV < 0.1"
x=t(as.matrix(groupAB_mean_diff_filt))
heatmap_image=heatmap.2(x, na.rm = TRUE, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=FALSE, dendrogram="row", main=main_title, ColSideColors=treatment_colors, cexRow=0.4, cexCol=0.8, col=bluered, symbreaks=TRUE)
dev.off()

#plot distributions of group mean values
pdf(file=file4)
par (mfrow=c(2,1))
hist(groupA_mean_data, col="blue", main="Synergy", xlab="Mean protein abundance", xlim=c(-2,2))
hist(groupB_mean_data, col="blue", main="Non-Synergy", xlab="Mean protein abundance", xlim=c(-2,2))
dev.off()



