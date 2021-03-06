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
data=rawdata
lines_meta=lines_meta_raw

#2. HER2+ lines only (Amplified only)
data=rawdata[rawdata[,"HER2_status"]=="Amplified",]
lines_meta=lines_meta_raw[lines_meta_raw[,"HER2_status"]=="Amplified",]

#choose either (A) or (B)
#Define comparison groups
#A. Synergy vs non-synergy
groupA_lines=which(lines_meta[,"Lap_AKT_Synergy"]=="Y")
groupB_lines=which(lines_meta[,"Lap_AKT_Synergy"]=="N")

#B. mutant (PI3K or PTEN) vs wildtype
groupA_lines=which(lines_meta[,"PI3K_OR_PTEN"]=="mutant")
groupB_lines=which(lines_meta[,"PI3K_OR_PTEN"]=="wt")

##Specify output files
#Choose only one of the following file names (1A, 1B, 2A, or 2B, etc)
#1A. All lines, synergy vs non-synergy 
wilcox_outfile="synergy_vs_nonsynergy_meantime_all_wilcox.txt"
pdf(file="synergy_vs_nonsynergy_all.pdf")
#pdf(file="synergy_vs_nonsynergy_all_vignettes.pdf")

#1B. All lines, mutant vs wild-type 
wilcox_outfile="mutant_vs_wildtype_meantime_all_wilcox.txt"
pdf(file="mutant_vs_wildtype_all.pdf")

#2A. HER2+ lines only (Amplified only), synergy vs non-synergy
wilcox_outfile="synergy_vs_nonsynergy_meantime_HER2amp_wilcox.txt"
pdf(file="synergy_vs_nonsynergy_HER2amp.pdf")
pdf(file="synergy_vs_nonsynergy_HER2amp_vignettes.pdf")

#2B. HER2+ lines only (Amplified only), mutant vs wild-type
wilcox_outfile="mutant_vs_wildtype_meantime_HER2amp_wilcox.txt"
pdf(file="mutant_vs_wildtype_HER2amp.pdf")

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

#Compile list of treatment/line combinations
treatment_line=vector(length=length(treatments)*length(lines))
n=0
for (i in 1:length(treatments)){
 for (j in 1:length(lines)){
  n=n+1
  treatment_line[n]=paste(lines[j],treatments[i],sep="_")
 }
}

treatment_line_array = array(0, dimnames = list(treatment_line, c("line","treatment")), dim=c(length(treatment_line),2))

#Mean protein abundance across all time points for each line/treatment
mean_data = array(0, dimnames = list(treatment_line, colnames(data)), dim=c(length(treatment_line),length(colnames(data))))
n=0
for (i in 1:length(treatments)){
 for (j in 1:length(lines)){
  n=n+1
  line_treatment_data=data[descriptor[,"line"]==lines[j] & descriptor[,"treatment"]==treatments[i],]
  mean_line_treatment_data=mean(line_treatment_data)
  mean_data[n,]=mean_line_treatment_data
  treatment_line_array[n,]=c(lines[j],treatments[i])
 }
}

#Now, perform comparisons between different subsets of cell lines for each treatment type for all proteins
proteins=colnames(mean_data)

#Wilcox comparison of cell line groupA vs groupB with abundances averaged across all time points
wilcox_results = array(0, dimnames = list(treatments, proteins), dim=c(length(treatments),length(proteins)))

for (i in 1:length(treatments)){
 treatment_mean_data=mean_data[treatment_line_array[,"treatment"]==treatments[i],]
 groupA_lines_data=treatment_mean_data[groupA_lines,]
 groupB_lines_data=treatment_mean_data[groupB_lines,]
 for (j in 1:length(proteins)){
  wilcox_result=wilcox.test(x=groupA_lines_data[,j], y=groupB_lines_data[,j], alternative="two.sided", paired=FALSE)
  wilcox_results[treatments[i],proteins[j]]=wilcox_result$p.value
 }
}

#Need to do multiple testing correction

#Write results to file
wilcox_results_out=t(rbind(protein_meta[1,],wilcox_results))
write.table (wilcox_results_out, sep="\t", file=wilcox_outfile)


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

#Create plots for each protein comparing mean abundance between groupA and groupB cell lines for each treatment
colors=rainbow(4)

#Choose one legend depending on group being analysed
#A. Synergy vs non-synergy
legend_text=c("Syn - DMSO","Syn - Lap","Syn - AKTi","Syn - Lap/AKTi","NonSyn - DMSO","NonSyn - Lap","NonSyn - AKTi","NonSyn - Lap/AKTi")

#B. mutant (PI3K or PTEN) vs wildtype
#legend_text=c("mt - DMSO","mt - Lap","mt - AKTi","mt - Lap/AKTi","wt - DMSO","wt - Lap","wt - AKTi","wt - Lap/AKTi")

ymin=min(c(groupA_mean_data,groupB_mean_data))
ymax=max(c(groupA_mean_data,groupB_mean_data))
xval=c(1:8)

par(mfrow=c(3,3))
#par(mfrow=c(1,1))

#Plot legend in first empty plot for reference
plot(x=c(1,8), y=c(ymin,ymax), type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
legend("topleft", legend=legend_text, lty=c(2,2,2,2,1,1,1,1), lwd=2, col=c(colors,colors), cex=0.9)

for (i in 1:length(proteins)){
 protein=proteins[i]
 plot(x=c(1,8), y=c(ymin,ymax), type="n", xaxt="n", xlab="time (hr)", ylab="protein abundance", main=paste(protein, protein_names[protein], sep=" - "), cex.main=0.75)
 axis(side=1, at=xval, labels=timepoints)

 lines(x=xval, y=groupA_mean_data[treatment_time_array[,"treatment"]=="DMSO",protein], col=colors[1], lty=2, lwd=2)
 lines(x=xval, y=groupA_mean_data[treatment_time_array[,"treatment"]=="Lapatinib",protein], col=colors[2], lty=2, lwd=2)
 lines(x=xval, y=groupA_mean_data[treatment_time_array[,"treatment"]=="AKTi",protein], col=colors[3], lty=2, lwd=2)
 lines(x=xval, y=groupA_mean_data[treatment_time_array[,"treatment"]=="Lap_AKTi",protein], col=colors[4], lty=2, lwd=2)

 lines(x=xval, y=groupB_mean_data[treatment_time_array[,"treatment"]=="DMSO",protein], col=colors[1], lwd=2)
 lines(x=xval, y=groupB_mean_data[treatment_time_array[,"treatment"]=="Lapatinib",protein], col=colors[2], lwd=2)
 lines(x=xval, y=groupB_mean_data[treatment_time_array[,"treatment"]=="AKTi",protein], col=colors[3], lwd=2)
 lines(x=xval, y=groupB_mean_data[treatment_time_array[,"treatment"]=="Lap_AKTi",protein], col=colors[4], lwd=2)

 text(1, -1.75, paste("RPPA Qual =",RPPA_qualities[protein], sep=" "), adj = c(0,0), col=RPPA_qual_col[protein][,])
}

dev.off()

