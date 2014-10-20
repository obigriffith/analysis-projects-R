library(multtest)

#Gene
genedatafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Matrix_GeneExpression_v53.txt"
geneexpressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Expressed_GeneExpression_v53.txt"

#Transcript
featuredatafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Matrix_TranscriptExpression_v53.txt"
featureexpressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Expressed_TranscriptExpression_v53.txt"

#KnownJunction
#featuredatafile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Matrix_KnownJunctionExpression_v53_test.txt"
#featureexpressedfile="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/matrix/Expressed_KnownJunctionExpression_v53_test.txt"

outdir="C:/Users/Obi/Documents/Projects/BCCL/ALEXA/"
outfile="ANCOVA_test.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_gene_data=read.table(genedatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
gene_exp_status=read.table(geneexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_feat_data=read.table(featuredatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
feat_exp_status=read.table(featureexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

setwd(outdir)
header=colnames(raw_gene_data)

pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 0.7 #Minimum coefficient of variation (0.7 recommended?)
cov_max = 10 #Maximum cov

groupMemberLists="Basal:HCC1806,HCC1143,SUM149PT,HCC1937,HCC1187,HCC1954,HCC3153,HCC1569,HCC70,HCC1500;Luminal:X600MPE,ZR75B,HCC1428,HCC202,MCF7,CAMA1,ZR7530,BT474,HCC1419,MDAMB13v1"
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
libs=vector(length=libcount)
lib_group=vector(length=libcount)
k=0
for (i in 1:length(group_libs)){
  for (j in 1:length(group_libs[[i]])){
    k=k+1
    libs[k]=group_libs[[i]][j]
    lib_group[k]=group_names[i]
  }
}

#Define a percent expressed function
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}

#Apply PE fun to both feature-level and gene-level
w=feat_exp_status[,libs]
pe_data=apply(w, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
feat_data=raw_feat_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}

#Apply COV function to feature expression only (cases where gene is invariant can still be interesting)
y=feat_data[,libs]
cov_data=apply(y, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
feat_data=feat_data[passed_cov,]
feat_data_filt_genes=feat_data[,"EnsEMBL_Gene_ID"]

gene_ids=raw_gene_data[,"EnsEMBL_Gene_ID"]
gene_data=raw_gene_data
rownames(gene_data)=gene_ids
gene_data=gene_data[feat_data_filt_genes,]

#Create dataframe starting with feature data
data=feat_data

#Grab just library data for genes/features remaining after filtering
z=feat_data[,libs] #feature data
z2=gene_data[,libs] #gene data

#Define function to calculate fold change - Fold changes will always be >= 1 (no change) and +ve versus -ve indicates the direction of the change
#Add '1' to all values to stabilize variance and prevent divion by 0
fold.change=function(x){
  dataA = x[group_libs[[1]]]
  dataB = x[group_libs[[2]]]
  mean_dataA = mean(dataA)
  mean_dataB = mean(dataB)

  #calculate fold change value
  if (mean_dataA >= mean_dataB){
      result = (mean_dataA+1)/(mean_dataB+1)
    }else{
      result = ((mean_dataB+1)/(mean_dataA+1))*-1
    }
  return(result)
}

#Define a function to calculate log2 difference
log2.diff=function(x){
  dataA = x[group_libs[[1]]]
  dataB = x[group_libs[[2]]]
  mean_dataA = mean(dataA)
  mean_dataB = mean(dataB)
  result = log2(mean_dataA+1) - log2(mean_dataB+1)
  return(result)
}

#Define a function to calculate SI value
SI_fun=function(x){
  A_SEQ_Norm = mean(as.numeric(z[x,group_libs[[1]]]))
  A_GENE_Norm = mean(as.numeric(z2[x,group_libs[[1]]]))
  B_SEQ_Norm = mean(as.numeric(z[x,group_libs[[2]]]))
  B_GENE_Norm = mean(as.numeric(z2[x,group_libs[[2]]]))
  SI = log2( ((A_SEQ_Norm+1)/(A_GENE_Norm+1)) / ((B_SEQ_Norm+1)/(B_GENE_Norm+1)) )
  return(SI)
}

#Define function to calculate ANOVA statistic
ANCOVA_fun=function(x){
  xdf=data.frame(t(z[x,]), t(z2[x,]), lib_group)
  names(xdf)=c("feat_expression", "gene_expression", "lib_group")
#  ANCOVA_result=aov(feat_expression ~ lib_group + gene_expression, data=xdf)
  ANCOVA_result=aov(feat_expression ~ gene_expression + lib_group, data=xdf)
#  ANCOVA_result1=aov(feat_expression ~ gene_expression + lib_group, data=xdf)
#  ANCOVA_result2=aov(feat_expression ~ gene_expression + lib_group + gene_expression:lib_group, data=xdf)
#  ANCOVA_result2=aov(feat_expression ~ gene_expression*lib_group, data=xdf)
  result=summary(ANCOVA_result)[[1]]["lib_group","Pr(>F)"]
#  lib_group_main=summary(ANCOVA_result1)[[1]]["lib_group","Pr(>F)"]
#  lib_group_main2=summary(ANCOVA_result2)[[1]]["lib_group","Pr(>F)"]
#  interaction=summary(ANCOVA_result2)[[1]]["gene_expression:lib_group","Pr(>F)"]
#  return(c(lib_group_main,lib_group_main2,interaction))
  return(result)
}

#Grab indexes of all remaining data (use to pull out corresponding feature/gene from two datafiles)
#index=as.matrix(rownames(z))
index=as.matrix(1:length(z[,1]))

#If number of groups is 2 calculate FC and log2 Diff on mean values, otherwise set to NA:
#Calculate fold change and log2 difference for both GENES and SEQS
#Then calculate SI value
if (length(groups)==2){
  print ("Calculating fold change values and log2 differences")
  data[,"SEQ_Fold_Change"] = apply(z, 1, fold.change)
  data[,"SEQ_Log2_Diff"] = apply(z, 1, log2.diff)
  data[,"SEQ_groupA_mean"] = apply(z[,group_libs[[1]]], 1, mean)
  data[,"SEQ_groupB_mean"] = apply(z[,group_libs[[2]]], 1, mean)
  data[,"GENE_Fold_Change"] = apply(z2, 1, fold.change)
  data[,"GENE_Log2_Diff"] = apply(z2, 1, log2.diff)
  data[,"GENE_groupA_mean"] = apply(z2[,group_libs[[1]]], 1, mean)
  data[,"GENE_groupB_mean"] = apply(z2[,group_libs[[2]]], 1, mean)
  data[,"SI"] = apply(index, 1, SI_fun)

  #Determine which events exhibit a reciprocal shift ((Gene DE +ve AND Seq DE -ve) OR (Gene DE -ve AND Seq DE +ve))
  data[,"Reciprocal"] = 0
  i = which((data[,"GENE_Log2_Diff"] > 0 & data[,"SEQ_Log2_Diff"] < 0) | (data[,"GENE_Log2_Diff"] < 0 & data[,"SEQ_Log2_Diff"] > 0))  
  data[i,"Reciprocal"] = 1

  #Calculate the degree of reciprocity
  #Reciprocity Score = abs((abs(GENE_DE) + abs(SEQ_DE)) / (abs(GENE_DE) -abs(SEQ_DE)))
  data[,"Reciprocity"] = (abs(data[,"GENE_Log2_Diff"]) + abs(data[,"SEQ_Log2_Diff"])) / (abs(data[,"GENE_Log2_Diff"]) - abs(data[,"SEQ_Log2_Diff"]))

  #Calculate the percentage of differential expression contributed by the SEQ relative to the GENE
  #percent_SEQ_DE = (abs(SEQ_DE)/(abs(GENE_DE)+abs(SEQ_DE)))*100
  data[,"percent_SEQ_Log2_DE"] = (abs(data[,"SEQ_Log2_Diff"]) / (abs(data[,"SEQ_Log2_Diff"]) + abs(data[,"GENE_Log2_Diff"])))*100

}else{
  data[,"SEQ_Fold_Change"] = "NA"
  data[,"SEQ_Log2_Diff"] = "NA"
  data[,"SEQ_groupA_mean"] = "NA"
  data[,"SEQ_groupB_mean"] = "NA"
  data[,"GENE_Fold_Change"] = "NA"
  data[,"GENE_Log2_Diff"] = "NA"
  data[,"GENE_groupA_mean"] = "NA"
  data[,"GENE_groupB_mean"] = "NA"
  data[,"SI"] = "NA"
  data[,"Reciprocal"] = "NA"
  data[,"Reciprocity"] = "NA"
  data[,"percent_SEQ_Log2_DE"] = "NA"
}

#Calculate ANCOVA p-values with feature expression as dependent variable and gene expression as a covariate
#Consider rank-transforming data to create non-parametric test
print ("Calculating ANCOVA p-values")

#Apply the ANCOVA function to all rows of the data for the appropriate columns
#ANCOVA_results=apply(index, 1, ANCOVA_fun)
ANCOVA_results=t(apply(index, 1, ANCOVA_fun))

#Correct p-values
print ("Correcting p-value for multiple testing")
ANCOVA_pvalues_adj=mt.rawp2adjp(as.numeric(ANCOVA_results), proc="BH")
ANCOVA_pvalues_adj_orig_order=ANCOVA_pvalues_adj$adjp[order(ANCOVA_pvalues_adj$index),]
data=cbind(data[,1:4],data[,c("SEQ_groupA_mean","SEQ_groupB_mean","SEQ_Fold_Change","SEQ_Log2_Diff","GENE_groupA_mean","GENE_groupB_mean","GENE_Fold_Change","GENE_Log2_Diff","SI","Reciprocal","Reciprocity","percent_SEQ_Log2_DE")], ANCOVA_pvalues_adj_orig_order)
#data=cbind(data[,1:4],data[,c("SEQ_groupA_mean","SEQ_groupB_mean","SEQ_Fold_Change","SEQ_Log2_Diff","GENE_groupA_mean","GENE_groupB_mean","GENE_Fold_Change","GENE_Log2_Diff","SI","Reciprocal","Reciprocity","percent_SEQ_Log2_DE")], ANCOVA_results)


#Write to file the log2 DE, foldchange, and SI values for every input element (gene, exon, junction, etc.) - also write a version sorted it on SI/pvalue
print ("Writing fold-change output text files")
write.table (data, sep="\t", file=outfile, quote=FALSE, row.names=FALSE)


