library(genefilter)
library(heatmap.plus)
library(fBasics)

#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop8/")

datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/processing/processed_final2/train_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)
clindatafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.train.2.txt"
combined_case_pred_file="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/Cepheid_CasePredictions_combined.txt"

top100heatmap_pdffile="top100heatmap_k8.pdf"
GeneClusterfile="gene_clusters.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")
case_predictions_all_combined_import=read.table(combined_case_pred_file, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
case_data_order=order(case_predictions_all_combined_import[,"GSM"])
case_predictions_all_combined=case_predictions_all_combined_import[case_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)

#If there are predictor variables that are constant/invariant, consider removing them
#Preliminary gene filtering
X=rawdata[,4:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=rawdata[filt,] 

#Get potential predictor variables
predictor_data=t(filt_Data[,4:length(header)]) #Filtered
predictor_names=paste(filt_Data[,3]," (",filt_Data[,1],")", sep="") #Filtered, gene symbol + probe ids
colnames(predictor_data)=predictor_names

#Print to screen and make sure samples line up correctly
cbind(rownames(predictor_data),as.vector(clindata[,"GSM"]),as.vector(case_predictions_all_combined[,"GSM"]))

#Filter down to just cases which have 10yr FU (i.e., exclude NAs)
cases_10yr = !is.na(clindata[,"X10yr_relapse"])
clindata_10yr=clindata[cases_10yr,]
predictor_data_10yr=predictor_data[cases_10yr,]
case_predictions_all_combined_10yr=case_predictions_all_combined[cases_10yr,]

#Get target variable and specify as factor/categorical
target=clindata_10yr[,"X10yr_relapse"] #recurrences after 10yrs not considered events
target[target==0]="NoRelapse"
target[target==1]="Relapse"
target=as.factor(target)

#Load top 100 genes/probes (pre-defined previously)
top100file="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/top100_Vars_nonRedundant.txt" #combined (standardCDF + customCDF)
#Read in data (expecting a tab-delimited file with header line and rownames)
top100_data_import=read.table(top100file, header = FALSE, na.strings = "NA", sep="\t")
rownames(top100_data_import)=top100_data_import[,4] #Set gene/probe names as rownames

#Exclude bad probes
bad_probes=c("KIAA0101 (9768_at)","PTTG1 (203554_x_at)","FRG1 (2483_at)","DPP3 (218567_x_at)","PDCD6 (222380_s_at)","KIAA0776 (212634_at)","TAF1D (218750_at)","SNORA25 (684959_at)","FOLR1 (211074_at)","KIAA1467 (57613_at)")
top100_data_import=top100_data_import[-which(rownames(top100_data_import) %in% bad_probes),]
top100=as.vector(top100_data_import[,4]) #gene/probe labels
top100_predictor_data_10yr=predictor_data_10yr[,top100]

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Perform k-means clustering to break into desired number of clusters
kval=8
#kmeans_clusters=kmeans(t(top100_predictor_data_10yr), kval)
kmeans_clusters=kmeans(t(top100_predictor_data_10yr), kval, nstart = 25)

#cluster kmeans clusters to get order to display the clusters
#kmeans_cluster_order=heatmap.2(as.matrix(kmeans_clusters$centers), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", labCol=NA, RowSideColors=rainbow(8), cexRow=1, col=rev(heat.colors(75)))$rowInd
kmeans_cluster_order=heatmap.plus(as.matrix(kmeans_clusters$centers), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", labCol=NA, cexRow=1, col=rev(heat.colors(75)))$rowInd

#Now get order of genes
genes_kmeans_order_new=rep(NA,length(kmeans_clusters$cluster))
topXgenes=rep(NA,length(kmeans_cluster_order))
for (k in 1:kval){
 cluster=kmeans_cluster_order[k]
 cluster_genes=kmeans_clusters$cluster[kmeans_clusters$cluster==cluster]
 i=length(which(!is.na(genes_kmeans_order_new)))+1
 j=i+length(cluster_genes)-1
 genes_kmeans_order_new[i:j]=names(cluster_genes)
 topXgenes[k]=names(cluster_genes[1])
}

#Create color side bar for heatmap based on kmeans clusters
colors=rainbow(kval)
kmeans_colors=colors[as.vector(kmeans_clusters$cluster)]
kmeans_colors_ordered=colors[as.vector(kmeans_clusters$cluster[genes_kmeans_order_new])]

top100_ordered=genes_kmeans_order_new #gene/probe labels - reordered according to kmeans

#Get actual variable importances, and create colors
top100VarImp=as.vector(top100_data_import[,3]) #ordered already from large to small
top100VarImp_colors=seqPalette(length(top100VarImp), name="Blues")[round(top100VarImp*length(top100VarImp))]
top100VarImp_ordered=as.vector(top100_data_import[genes_kmeans_order_new,3]) #kmeans order
top100VarImp_colors_ordered=seqPalette(100, name="Blues")[round(top100VarImp_ordered*100)]

#Create dataframe with all genes and their corresponding cluster and Var Imps
gene_cluster_data=cbind(genes_kmeans_order_new,as.vector(kmeans_clusters$cluster[genes_kmeans_order_new]),top100VarImp_ordered)
colnames(gene_cluster_data)=c("gene","cluster","VarImp")
write.table(gene_cluster_data, file=GeneClusterfile, sep="\t", row.names=FALSE)

#Print out k best genes in same order as they will appear on heatmap (from top to bottom)
getBestClusterGene=function(x){gene_cluster_data[which(gene_cluster_data[,2]==x),][1,1]}
sapply(rev(kmeans_cluster_order),getBestClusterGene)



#Set target colors for patients
target_colors=as.vector(target)
target_colors[target_colors=="NoRelapse"]="yellow"
target_colors[target_colors=="Relapse"]="red"

risk_group_colors=as.vector(case_predictions_all_combined_10yr[,"RF_Group1"])
risk_group_colors[risk_group_colors=="low"]="lightblue"
risk_group_colors[risk_group_colors=="int"]="blue"
risk_group_colors[risk_group_colors=="high"]="darkblue"

#Get actual relapse risk scores so that heatmap can be ordered with this
risk_group_scores=case_predictions_all_combined_10yr[,"Relapse"]
risk_order=order(risk_group_scores)
risk_group_colors_ordered=risk_group_colors[risk_order]
target_colors_ordered=target_colors[risk_order]


#Reformat with other heatmap function to allow multiple color side bars
clab=cbind(target_colors,risk_group_colors)
clab_riskorder=cbind(target_colors_ordered,risk_group_colors_ordered)
colnames(clab)=c("Relapse status","Risk group")
colnames(clab_riskorder)=c("Relapse status","Risk group")
rlab=cbind(kmeans_colors, top100VarImp_colors)
rlab_kmeansorder=cbind(kmeans_colors_ordered, top100VarImp_colors_ordered)
colnames(rlab)=c("Kmeans cluster","Var. Imp.")
colnames(rlab_kmeansorder)=c("Kmeans cluster","Var. Imp.")

pdf(file=top100heatmap_pdffile)
#Create separate heatmaps. One with columns/rows ordered by hclust, the other by risk score and kmeans cluster
#Set par margins for heatmap.plus (heatmap.2 has to be handled differently)?
par(oma=c(3,4.5,1,2)) #bottom, left, top, right

#Will need to rearrange the actual patient/gene data and set Colv/Rowv=FALSE because reordering with a vector doesn't seem to work as expected.
top100_predictor_data_10yr_riskorder=top100_predictor_data_10yr[risk_order,genes_kmeans_order_new]
heatmap.plus(as.matrix(t(top100_predictor_data_10yr_riskorder)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Colv=NA, Rowv=NA, labRow=top100_ordered, labCol=NA, ColSideColors=clab_riskorder, RowSideColors=rlab_kmeansorder, cexRow=0.50, col=rev(heat.colors(75)))
dev.off()


