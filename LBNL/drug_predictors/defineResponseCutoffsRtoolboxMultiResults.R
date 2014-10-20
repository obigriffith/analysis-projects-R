library("ggplot2")
library("mclust")

#Load probability matrix for drugs
probmatrixfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/Results/TCGA_results_ExpMethCNV_prob_matrix.txt"
probmatrix=read.table(probmatrixfile,sep="\t",header=TRUE, row.names=1)
probability_cutoffs_file="TCGA_results_ExpMethCNV_prob_cutoffs.txt"
probability_classes_file="TCGA_results_ExpMethCNV_prob_classes.txt"
probability_cutoffs_plots_file="TCGA_results_ExpMethCNV_prob_cutoff_plots.pdf"
probability_classes_filtered_file="TCGA_results_ExpMethCNV_prob_classes_filtered.txt"

#Set output dir and files
output_dir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/TCGA/RtoolboxResults/Results/"
setwd(output_dir)

drugs=colnames(probmatrix)
patients=rownames(probmatrix)

#Create matrix to store cutoffs and numbers of patients in each response group
cutoffs=matrix(NA, nrow=length(drugs), ncol=5, dimnames=list(drugs,c("cutoff1","cutoff2","LowN","IntN","HighN")))

#Create matrix to store probability risk classes
probclassmatrix=probmatrix

#For each drug, choose cutoffs and plot distributions
pdf(file=probability_cutoffs_plots_file, height=7.5, width=10)
#par(mfrow=c(2,2))
#layout(matrix(c(1,1,0,2), 2, 2, byrow=TRUE), respect=TRUE)
for (i in 1:length(drugs)){
#for (i in 1:6){
  #i=2
  drug=drugs[i]
  #Attempt to separate probabilities based on gaussian mixture model
  #Allow mclust to determine best model
  x=as.numeric(probmatrix[,drug])
  mclust_prob=Mclust(x)
  #summary(mclust_prob, x) #Gives you list of values returned
  classification_prob=mclust_prob$classification
  num_clusters=mclust_prob$G

  #Choose cutoffs - use max and min of left and right clusters to break into components
  if (num_clusters==1){ #If only a single distribution detected, use tertiles as cutoffs
#    tertiles=quantile(x, probs=c(1/3,2/3))
    tertiles=quantile(x, probs=c(1/4,3/4)) #Instead divide into bottom quartile, middle two quartiles, and top quartile
    cutoff1=tertiles[1]
    cutoff2=tertiles[2]
  }
  if (num_clusters==2){ #If two clusters just use one cutoff between them and median of larger cluster as second cutoff
    cluster_sizes=c(length(x[classification_prob==1]),length(x[classification_prob==2]))
    larger_cluster=which(cluster_sizes==max(cluster_sizes))
    if (larger_cluster==1){
      cutoff1=median(x[classification_prob==1])
      cutoff2=max(x[classification_prob==1])
    }else{
      cutoff1=min(x[classification_prob==2])
      cutoff2=median(x[classification_prob==2])
    }
    #If component #1 is within component #2 (i.e., not broken into left and right, but middle and outsides) then use that middle component to set cutoffs
    if (max(x[classification_prob==1]) < max(x[classification_prob==num_clusters]) & min(x[classification_prob==1]) > min(x[classification_prob==num_clusters])){
      cutoff1=min(x[classification_prob==1])
      cutoff2=max(x[classification_prob==1])
    }
  }
  if (num_clusters==3){ #If three clusters use max of 1st cluster and min of 3rd cluster to split into low/int/high groups
    cutoff1=max(x[classification_prob==1])
    cutoff2=min(x[classification_prob==num_clusters])
  }
  if (num_clusters==4){ #If four clusters use max of 2nd cluster and min of last cluster to split into low/int/high groups
    cutoff1=max(x[classification_prob==2])
    cutoff2=min(x[classification_prob==num_clusters])
  }
  if (num_clusters==5){ #If five clusters use max of 1st cluster and min of 2nd last cluster to split into low/int/high groups
    cutoff1=max(x[classification_prob==1])
    cutoff2=min(x[classification_prob==(num_clusters-1)])
  }
  if (num_clusters==6){ #If six clusters use max of 3rd cluster and min of 2nd last cluster to split into low/int/high groups
    cutoff1=max(x[classification_prob==3])
    cutoff2=min(x[classification_prob==(num_clusters-1)])
  }
  if (num_clusters==7){ #If seven clusters use max of 3rd cluster and min of 2nd last cluster to split into low/int/high groups
    cutoff1=max(x[classification_prob==4])
    cutoff2=min(x[classification_prob==(num_clusters-1)])
  }

  if (max(x)<0.65){ #If no patients have probability > 0.65, set all to int category
    cutoff1=0
    cutoff2=1
  }


  #Plot distribution of values for clusters and cutoff values
  data=data.frame(x,as.factor(classification_prob))
  colnames(data)=c("prob","cluster")

  if (cutoff1==0 & cutoff2==1){ #Don't plot cutoff lines
    p=ggplot(data, aes(prob, fill=cluster)) + geom_density(alpha = 0.2, stat="bin", binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30) + xlab("Response Prob") + ylab("Density") + opts(title=drug) + opts(legend.position = "none")
  }
  if (cutoff1==0 & cutoff2!=1){
    p=ggplot(data, aes(prob, fill=cluster)) + geom_density(alpha = 0.2, stat="bin", binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30) + xlab("Response Prob") + ylab("Density") + opts(title=drug) + geom_vline(xintercept = cutoff2) + opts(legend.position = "none")
  }
  if (cutoff1!=0 & cutoff2==1){
    p=ggplot(data, aes(prob, fill=cluster)) + geom_density(alpha = 0.2, stat="bin", binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30) + xlab("Response Prob") + ylab("Density") + opts(title=drug) + geom_vline(xintercept = cutoff1) + opts(legend.position = "none")
  }
  if (cutoff1!=0 & cutoff2!=1){
    p=ggplot(data, aes(prob, fill=cluster)) + geom_density(alpha = 0.2, stat="bin", binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30) + xlab("Response Prob") + ylab("Density") + opts(title=drug) + geom_vline(xintercept = cutoff1) + geom_vline(xintercept = cutoff2)+ opts(legend.position = "none")
  }

  #p=ggplot(data, aes(prob, colour=cluster)) + geom_freqpoly(binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30) + xlab("Response Prob") + ylab("Density") + opts(title=drug) + geom_vline(xintercept = cutoff1) + geom_vline(xintercept = cutoff2)
  #p=ggplot(data, aes(prob)) + geom_area(aes(y = ..count.., fill = cluster, group = cluster), binwidth=(max(data[,"prob"])-min(data[,"prob"]))/30, alpha=0.3, position = "identity", stat = "bin") + opts(axis.text.x=theme_text(angle=-60, hjust=0, size = 6), legend.position = "none") + geom_vline(xintercept = cutoff1) + geom_vline(xintercept = cutoff2)
  print(p)

  prob_low=which(x <= cutoff1)
  prob_int=which(x > cutoff1 & x < cutoff2)
  prob_high=which(x >= cutoff2)

  prob_low_count=length(prob_low)
  prob_int_count=length(prob_int)
  prob_high_count=length(prob_high)

  prob_classes=x
  prob_classes[prob_low]="res"
  prob_classes[prob_int]="int"
  prob_classes[prob_high]="sens"

  cutoffs[drug,"cutoff1"]=cutoff1
  cutoffs[drug,"cutoff2"]=cutoff2
  cutoffs[drug,"LowN"]=prob_low_count
  cutoffs[drug,"IntN"]=prob_int_count
  cutoffs[drug,"HighN"]=prob_high_count
  probclassmatrix[,drug]=prob_classes

}
dev.off()
write.table(cbind(drugs,cutoffs), file=probability_cutoffs_file, row.names=FALSE, sep="\t")
write.table(cbind(patients,probclassmatrix), file=probability_classes_file, row.names=FALSE, sep="\t")

#Create subset file for drugs passing filter
filtered_drugs=c("X5.FU","AG1478","AKT1.2.inhibitor","BIBW2992","Bortezomib.batch1","Cisplatin","Etoposide","Fascaplysin","GSK.MEKi","GSK2..PLKi.","GSK2119563A","GSK2126458A","GSK2141795c","GSK615B.PI3Ki.","GSK_Tykerb","Iressa","Ixabepilone","NU6102","Nutlin.3a","PF.3084014","Paclitaxel","Tamoxifen")
probclassmatrix_filtered=probclassmatrix[,filtered_drugs]
write.table(cbind(patients,probclassmatrix_filtered), file=probability_classes_filtered_file, row.names=FALSE, sep="\t")

