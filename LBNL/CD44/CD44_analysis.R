library(multtest)
library(grid)

#Expression data
datafile="C:/Users/Obi/Documents/Projects/CD44/ALEXA/25lines/ENSG00000026508.txt"

outdir="C:/Users/Obi/Documents/Projects/CD44/ALEXA/25lines/"
outfile="CD44_analysis.txt"
pdffile="CD44_analysis.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

setwd(outdir)
header=colnames(raw_data)
features=raw_data[,2]

#Set feature names as rowname
data=raw_data
rownames(data)=features

#Specify subtype for each library
groupMemberLists="Basal:HCC1806,HCC1143,SUM149PT,HCC1937,HCC1954,HCC3153,HCC1569,HCC70,HCC1500;Luminal:X600MPE,ZR75B,HCC1428,HCC202,MCF7,CAMA1,ZR7530,BT474,HCC1419,MDAMB134VI;Basal_NM:HCC1395,MDAMB231;ClaudinLow:MCF10A,MCF10F;Normal:Normal;Other:M4A4;Failed:HCC1187"
groups=strsplit(groupMemberLists, split=";")[[1]]
group_libs=vector("list", length=length(groups))
group_names=vector(length=length(groups))

libcount=0
for (i in 1:length(groups)){
  group_name=strsplit(groups[i],split=":")[[1]][1]
  libraries_str=strsplit(groups[i],split=":")[[1]][2]
  libraries=strsplit(libraries_str,split=",")[[1]]
  libcount=libcount+length(libraries)
  group_names[i]=group_name
  group_libs[[i]]=libraries
}

#Make vectors of libraries and their groups in the same order for look up purposes
libs_all=vector(length=libcount)
lib_group_all=vector(length=libcount)
k=0
for (i in 1:length(group_libs)){
  for (j in 1:length(group_libs[[i]])){
    k=k+1
    libs_all[k]=group_libs[[i]][j]
    lib_group_all[k]=group_names[i]
  }
}

#Define function to calculate ANOVA statistic for DE
ANOVA_fun=function(x){
  xdf=data.frame(x,lib_group)
  names(xdf)=c("expression", "lib_group")
  ANOVA_result=aov(expression~lib_group, data=xdf)
  result=summary(ANOVA_result)[[1]]["lib_group","Pr(>F)"]
  return(result)
}

#Define function to calculate ANCOVA statistic for AE
ANCOVA_fun=function(x){
  xdf=data.frame(x, t(z2), lib_group)
  names(xdf)=c("feat_expression", "gene_expression", "lib_group")
  ANCOVA_result=aov(feat_expression~gene_expression*lib_group, data=xdf)
  result=summary(ANCOVA_result)[[1]]["lib_group","Pr(>F)"]
  return(result)
}

#Define a function for doing MTC
MTC_fun=function(x){
  pvalues_adj=mt.rawp2adjp(as.numeric(x), proc="BH")
  pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
  return(pvalues_adj_orig_order)
}



#Select subset for analysis (Choose only one from below)
#1. All vs All
#comp_name="All_vs_All"
#subset=1:length(lib_group_all)

#2. Basal vs Luminal
#comp_name="Basal_vs_Luminal"
#subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal")

#3. All vs ALL (excluding failed)
#comp_name="All_vs_All_nofail"
#subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow" | lib_group_all=="Normal" | lib_group_all=="Other")

#4. Basal vs Basal_NM vs Luminal vs ClaudinLow
comp_name="Basal_vs_BasalNM_vs_Luminal_vs_ClaudinLow"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]

#Grab just library data for gene and all features
z=data[,libs]+1 #gene + feature data
z2=data[1,libs]+1 #gene data only

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)
ANOVA_results_adj=MTC_fun(ANOVA_results)
colnames(ANOVA_results_adj)=c(paste(comp_name,"AN","rawp",sep="_"),paste(comp_name,"AN","BH",sep="_"))

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(z, 1, ANCOVA_fun)
ANCOVA_results_adj=MTC_fun(ANCOVA_results)
colnames(ANCOVA_results_adj)=c(paste(comp_name,"ANC","rawp",sep="_"),paste(comp_name,"ANC","BH",sep="_"))

data=cbind(data, ANOVA_results_adj, ANCOVA_results_adj)


#Or, choose a subset from below and re-organize group names for targetted comparisons
#5. Basal vs (Basal_NM, Luminal, ClaudinLow)
comp_name="Basal_vs_Other"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]
lib_group[lib_group=="Luminal" | lib_group=="Basal_NM" | lib_group=="ClaudinLow"]="Other"

#Grab just library data for gene and all features
z=data[,libs]+1 #gene + feature data
z2=data[1,libs]+1 #gene data only

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)
ANOVA_results_adj=MTC_fun(ANOVA_results)
colnames(ANOVA_results_adj)=c(paste(comp_name,"AN","rawp",sep="_"),paste(comp_name,"AN","BH",sep="_"))

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(z, 1, ANCOVA_fun)
ANCOVA_results_adj=MTC_fun(ANCOVA_results)
colnames(ANCOVA_results_adj)=c(paste(comp_name,"ANC","rawp",sep="_"),paste(comp_name,"ANC","BH",sep="_"))

#Calculate FC values
FC = apply(z[,lib_group=="Basal"],1,mean) / apply(z[,lib_group=="Other"],1,mean)
FC[FC<1 & !is.na(FC)]=-1/FC[FC<1 & !is.na(FC)]

data=cbind(data, ANOVA_results_adj, ANCOVA_results_adj, FC)


#6. Basal_NM vs (Basal, Luminal, ClaudinLow)
comp_name="BasalNM_vs_Other"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]
lib_group[lib_group=="Basal" | lib_group=="Luminal" | lib_group=="ClaudinLow"]="Other"

#Grab just library data for gene and all features
z=data[,libs]+1 #gene + feature data
z2=data[1,libs]+1 #gene data only

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)
ANOVA_results_adj=MTC_fun(ANOVA_results)
colnames(ANOVA_results_adj)=c(paste(comp_name,"AN","rawp",sep="_"),paste(comp_name,"AN","BH",sep="_"))

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(z, 1, ANCOVA_fun)
ANCOVA_results_adj=MTC_fun(ANCOVA_results)
colnames(ANCOVA_results_adj)=c(paste(comp_name,"ANC","rawp",sep="_"),paste(comp_name,"ANC","BH",sep="_"))

#Calculate FC values
FC = apply(z[,lib_group=="Basal_NM"],1,mean) / apply(z[,lib_group=="Other"],1,mean)
FC[FC<1 & !is.na(FC)]=-1/FC[FC<1 & !is.na(FC)]

data=cbind(data, ANOVA_results_adj, ANCOVA_results_adj, FC)


#7. Luminal vs (Basal, Basal_NM, ClaudinLow)
comp_name="Luminal_vs_Other"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]
lib_group[lib_group=="Basal" | lib_group=="Basal_NM" | lib_group=="ClaudinLow"]="Other"

#Grab just library data for gene and all features
z=data[,libs]+1 #gene + feature data
z2=data[1,libs]+1 #gene data only

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)
ANOVA_results_adj=MTC_fun(ANOVA_results)
colnames(ANOVA_results_adj)=c(paste(comp_name,"AN","rawp",sep="_"),paste(comp_name,"AN","BH",sep="_"))

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(z, 1, ANCOVA_fun)
ANCOVA_results_adj=MTC_fun(ANCOVA_results)
colnames(ANCOVA_results_adj)=c(paste(comp_name,"ANC","rawp",sep="_"),paste(comp_name,"ANC","BH",sep="_"))

#Calculate FC values
FC = apply(z[,lib_group=="Luminal"],1,mean) / apply(z[,lib_group=="Other"],1,mean)
FC[FC<1 & !is.na(FC)]=-1/FC[FC<1 & !is.na(FC)]

data=cbind(data, ANOVA_results_adj, ANCOVA_results_adj, FC)


#8. ClaudinLow vs (Basal, Basal_NM, Luminal)
comp_name="ClaudinLow_vs_Other"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]
lib_group[lib_group=="Basal" | lib_group=="Basal_NM" | lib_group=="Luminal"]="Other"

#Grab just library data for gene and all features
z=data[,libs]+1 #gene + feature data
z2=data[1,libs]+1 #gene data only

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)
ANOVA_results_adj=MTC_fun(ANOVA_results)
colnames(ANOVA_results_adj)=c(paste(comp_name,"AN","rawp",sep="_"),paste(comp_name,"AN","BH",sep="_"))

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(z, 1, ANCOVA_fun)
ANCOVA_results_adj=MTC_fun(ANCOVA_results)
colnames(ANCOVA_results_adj)=c(paste(comp_name,"ANC","rawp",sep="_"),paste(comp_name,"ANC","BH",sep="_"))

#Calculate FC values
FC = apply(z[,lib_group=="ClaudinLow"],1,mean) / apply(z[,lib_group=="Other"],1,mean)
FC[FC<1 & !is.na(FC)]=-1/FC[FC<1 & !is.na(FC)]

data=cbind(data, ANOVA_results_adj, ANCOVA_results_adj, FC)


#Write results to file
write.table (data, sep="\t", file=outfile, quote=FALSE, row.names=FALSE)

#Create boxplots comparing gene/feature expression levels for different subtypes
comp_name="Basal_vs_BasalNM_vs_Luminal_vs_ClaudinLow"
subset=which(lib_group_all=="Basal" | lib_group_all=="Luminal" | lib_group_all=="Basal_NM" | lib_group_all=="ClaudinLow")
libs=libs_all[subset]
lib_group=lib_group_all[subset]

#Grab just library data for gene and all features
z=t(data[,libs]) #gene + feature data
pdf(file=pdffile)
par(mfrow=c(3,3), las=2)
for (i in 1:length(features)){
  group_labels=c("Basal","ClaudinLow","BasalNM","Luminal")
  pvalue=formatC(data[features[i],"Basal_vs_BasalNM_vs_Luminal_vs_ClaudinLow_AN_BH"], digits=4, format="f")

  if(!is.na(data[features[i],"Basal_vs_Other_AN_BH"])){if(data[features[i],"Basal_vs_Other_AN_BH"]<0.05){group_labels[1]=paste(group_labels[1],"*",sep="")}}
  if(!is.na(data[features[i],"ClaudinLow_vs_Other_AN_BH"])){if(data[features[i],"ClaudinLow_vs_Other_AN_BH"]<0.05){group_labels[2]=paste(group_labels[2],"*",sep="")}}
  if(!is.na(data[features[i],"BasalNM_vs_Other_AN_BH"])){if(data[features[i],"BasalNM_vs_Other_AN_BH"]<0.05){group_labels[3]=paste(group_labels[3],"*",sep="")}}
  if(!is.na(data[features[i],"Luminal_vs_Other_AN_BH"])){if(data[features[i],"Luminal_vs_Other_AN_BH"]<0.05){group_labels[4]=paste(group_labels[4],"*",sep="")}}

  if(!is.na(data[features[i],"Basal_vs_Other_AN_BH"])){if(data[features[i],"Basal_vs_Other_ANC_BH"]<0.05){group_labels[1]=paste(group_labels[1],"+",sep="")}}
  if(!is.na(data[features[i],"ClaudinLow_vs_Other_AN_BH"])){if(data[features[i],"ClaudinLow_vs_Other_ANC_BH"]<0.05){group_labels[2]=paste(group_labels[2],"+",sep="")}}
  if(!is.na(data[features[i],"BasalNM_vs_Other_AN_BH"])){if(data[features[i],"BasalNM_vs_Other_ANC_BH"]<0.05){group_labels[3]=paste(group_labels[3],"+",sep="")}}
  if(!is.na(data[features[i],"Luminal_vs_Other_AN_BH"])){if(data[features[i],"Luminal_vs_Other_ANC_BH"]<0.05){group_labels[4]=paste(group_labels[4],"+",sep="")}}

  boxplot(z[lib_group=="Basal",i], z[lib_group=="ClaudinLow",i], z[lib_group=="Basal_NM",i], z[lib_group=="Luminal",i], names=group_labels, main=features[i], col=rainbow(4), cex.axis=0.9)
  legend("topright", legend=pvalue, box.lty=0)
}
dev.off()





