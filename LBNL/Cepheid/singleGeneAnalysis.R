library(mclust)
library("heatmap.plus")
library(gplots)
library(survival)

datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_train_survival/combined/ALL_gcrma.txt"
#datafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/processing/processed_final2/test_train_survival/customCDF/ALL_gcrma.txt"
clindatafile="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.txt"

#gene="LRIG1"
#gene="LARGE"
#gene="NCOA3" #AKA: AIB1
#gene="ESR1"
#gene="PADI2"
#gene="PADI4"
#gene="STAT1"
gene="PRLR"

outdir=paste("/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/single_gene/",gene,sep="")
outfile="test_train_survival_ALL_gcrma_gene_w_clin_data.txt"
outfile2="test_train_survival_ALL_gcrma_gene_mix_model.pdf"
outfile3="test_train_survival_ALL_gcrma_gene_hist.pdf"
outfile3B="test_train_survival_ALL_gcrma_gene_hist_tert.pdf"
outfile4="test_train_survival_ALL_gcrma_gene_heatmap.pdf"
survival_outfile="test_train_survival_ALL_gcrma_gene_KM_MM1.pdf"
survival_outfile2="test_train_survival_ALL_gcrma_gene_KM_MM2.pdf"
survival_outfile3="test_train_survival_ALL_gcrma_gene_KM_MM3.pdf"
survival_outfile4="test_train_survival_ALL_gcrma_gene_KM_TERT.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clin_data=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
raw_data=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

#Change to output dir
setwd(outdir)
header=colnames(raw_data)

#Define function to get single-gene Data
getGeneData=function(x) {
 Gene_Data_all=raw_data[which(raw_data[,3]==x),]
 Gene_exp_data=Gene_Data_all[,4:length(header)]
 probe_names=Gene_Data_all[,1]
 gene_names=Gene_Data_all[,3]
 names=paste(probe_names,"(",gene_names,")",sep="")
 rownames(Gene_exp_data)=names
 return(Gene_exp_data)
}

#Get gene data for gene of interest
gene_exp_data=getGeneData(gene)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap (single color sidebar)
x=as.matrix(gene_exp_data)
row_names=rownames(x)

pdf(file=outfile4)
heatmap.2(x, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="column", Rowv=FALSE, labRow=row_names, cexRow=0.5, labCol=FALSE, col=rev(heat.colors(75)))
dev.off()

#Determine best probe based on COV value
cov=apply(gene_exp_data,1,function(x){sd(x)/mean(x)})
best_probe=which(cov==max(cov))
gene_exp_data_best=gene_exp_data[best_probe,]
probe=rownames(gene_exp_data_best)

#Attempt to separate gene data based on gaussian mixture model
#Allow mclust to determine best model
x=as.numeric(gene_exp_data_best[1,])
mclust_gene=Mclust(x)
summary(mclust_gene, x) #Gives you list of values returned
classification_gene=mclust_gene$classification
num_clusters=mclust_gene$G

pdf(file=outfile2)
par(mfrow=c(num_clusters,1), oma=c(2,2,2,2))
colors=rainbow(num_clusters)

for (i in 1:num_clusters){
  hist(x[classification_gene==i], xlim=c(2,15), col=colors[i], xlab="Log2 GCRMA value", main=paste("component",i,sep=" "))
}
title(main="Separation of gene expression by model-based clustering", outer=TRUE)
dev.off()

#Choose cutoffs - use max and min of left and right clusters to break into components
cutoff1=max(x[classification_gene==1])
cutoff2=min(x[classification_gene==num_clusters])
#If component #1 is within component #2 (i.e., not broken into left and right, but middle and outsides) then use that middle component to set cutoffs
if (max(x[classification_gene==1]) < max(x[classification_gene==num_clusters]) & min(x[classification_gene==1]) > min(x[classification_gene==num_clusters])){
 cutoff1=min(x[classification_gene==1])
 cutoff2=max(x[classification_gene==1])
}

gene_low_count=length(which(x <= cutoff1))
gene_int_count=length(which(x > cutoff1 & x < cutoff2))
gene_high_count=length(which(x >= cutoff2))

#Also set cutoffs by tertile
tertiles=quantile(x, probs=c(1/3,2/3))
tert_cutoff1=tertiles[1]
tert_cutoff2=tertiles[2]
gene_tert1_count=length(which(x <= tert_cutoff1))
gene_tert2_count=length(which(x > tert_cutoff1 & x < tert_cutoff2))
gene_tert3_count=length(which(x >= tert_cutoff2))

#Histogram
#mixed-model
pdf(file=outfile3)
hist(as.numeric(x), xlim=c(2,15), breaks=40, col="blue", main=paste(probe, "histogram for all 858 samples", sep=" "), xlab="Log2 GCRMA value", ylab="Frequency")
abline(v=cutoff1, col="red")
abline(v=cutoff2, col="red")
legend("left", legend=c(paste("low cutoff=",cutoff1,sep=""),paste("high cutoff=",cutoff2,sep=""),paste("N_low=",gene_low_count,sep=""),paste("N_int=",gene_int_count,sep=""),paste("N_high=",gene_high_count,sep="")),bty="n")
dev.off()

#tertiles
pdf(file=outfile3B)
hist(as.numeric(x), xlim=c(2,15), breaks=40, col="blue", main=paste(probe, "histogram for all 858 samples", sep=" "), xlab="Log2 GCRMA value", ylab="Frequency")
abline(v=tert_cutoff1, col="red")
abline(v=tert_cutoff2, col="red")
legend("left", legend=c(paste("low cutoff=",format(tert_cutoff1, digits=3),sep=""),paste("high cutoff=",format(tert_cutoff2, digits=3),sep=""),paste("N_low=",gene_tert1_count,sep=""),paste("N_int=",gene_tert2_count,sep=""),paste("N_high=",gene_tert3_count,sep="")),bty="n")
dev.off()

#Extract low, int, high samples
gene_low=which(x <= cutoff1)
gene_int=which(x > cutoff1 & x < cutoff2)
gene_high=which(x >= cutoff2)

gene_tert1=which(x <= tert_cutoff1)
gene_tert2=which(x > tert_cutoff1 & x < tert_cutoff2)
gene_tert3=which(x >= tert_cutoff2)

#Join clinical and gene data and write to file
clin_gene_data=cbind(clin_data,t(gene_exp_data_best))

#Create new columns to store cutoff groups 
#Mixed-model groups
clin_gene_data[gene_low,"exp_group_mm"]="low"
clin_gene_data[gene_int,"exp_group_mm"]="int"
clin_gene_data[gene_high,"exp_group_mm"]="high"

#Mixed-model groups low/int vs high
clin_gene_data[gene_low,"exp_group_mm2"]="low"
clin_gene_data[gene_int,"exp_group_mm2"]="low"
clin_gene_data[gene_high,"exp_group_mm2"]="high"

#Mixed-model groups low vs int/high
clin_gene_data[gene_low,"exp_group_mm3"]="low"
clin_gene_data[gene_int,"exp_group_mm3"]="high"
clin_gene_data[gene_high,"exp_group_mm3"]="high"

#Tertiles
clin_gene_data[gene_tert1,"exp_group_tert"]="low"
clin_gene_data[gene_tert2,"exp_group_tert"]="int"
clin_gene_data[gene_tert3,"exp_group_tert"]="high"

#Create column for e_rfs where events beyond 10years are censored
clin_gene_data[,"e_rfs_10yrcens"]=clin_gene_data[,"e_rfs"]
clin_gene_data[clin_gene_data[,"t_rfs"] >= 10,"e_rfs_10yrcens"]=0

write.table(clin_gene_data, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#Association statistics
noRelapse_exp_values=clin_gene_data[clin_gene_data[,"X10yr_relapse"]=="0",probe]
Relapse_exp_values=clin_gene_data[clin_gene_data[,"X10yr_relapse"]=="1",probe]
wilcox_result=wilcox.test(x=noRelapse_exp_values, y=Relapse_exp_values, alternative="two.sided", paired=FALSE)


#Create survival plots
#Create new dataframe with just necessary data
surv_data=clin_gene_data[,c("t_rfs","e_rfs_10yrcens","exp_group_mm","exp_group_mm2","exp_group_mm3","exp_group_tert")]

#create a survival object using data
surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))


#Plot KM curve - Expression group Mixed-Model
#create kaplan-meier estimate, using gender as strata
krfit.by_exp_group_mm = survfit(surv_data.surv ~ exp_group_mm, data = surv_data)

#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ exp_group_mm, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

#plot K-M estimate
pdf(file=survival_outfile)
colors = rainbow(3)
title=paste("Survival by",gene,"expression (MM group)")
plot(krfit.by_exp_group_mm, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
groups=unique(surv_data[,"exp_group_mm"])
group_sizes=table(surv_data[,"exp_group_mm"])[groups]
legend_text=c(paste(groups, " ", "(", group_sizes, ")", sep=""),paste("p =", survpvalue, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()

#Plot KM curve - Expression group Mixed-Model 2
#create kaplan-meier estimate, using gender as strata
krfit.by_exp_group_mm2 = survfit(surv_data.surv ~ exp_group_mm2, data = surv_data)

#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ exp_group_mm2, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

#plot K-M estimate
pdf(file=survival_outfile2)
colors = rainbow(2)
title=paste("Survival by",gene,"expression (MM group: low/int vs high)")
plot(krfit.by_exp_group_mm2, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
groups=unique(surv_data[,"exp_group_mm2"])
group_sizes=table(surv_data[,"exp_group_mm2"])[groups]
legend_text=c(paste(groups, " ", "(", group_sizes, ")", sep=""),paste("p =", survpvalue, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()


#Plot KM curve - Expression group Mixed-Model 3
#create kaplan-meier estimate, using gender as strata
krfit.by_exp_group_mm3 = survfit(surv_data.surv ~ exp_group_mm3, data = surv_data)

#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ exp_group_mm3, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

#plot K-M estimate
pdf(file=survival_outfile3)
colors = rainbow(2)
title=paste("Survival by",gene,"expression (MM group: low vs int/high)")
plot(krfit.by_exp_group_mm3, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
groups=unique(surv_data[,"exp_group_mm3"])
group_sizes=table(surv_data[,"exp_group_mm3"])[groups]
legend_text=c(paste(groups, " ", "(", group_sizes, ")", sep=""),paste("p =", survpvalue, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()


#Plot KM curve - Expression group tertiles
#create kaplan-meier estimate, using gender as strata
krfit.by_exp_group_tert = survfit(surv_data.surv ~ exp_group_tert, data = surv_data)

#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ exp_group_tert, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = round(as.numeric(survpvalue), digits=3)

#Calculate p-value for linear test
surv_data_lin=clin_gene_data[,c("t_rfs","e_rfs_10yrcens","exp_group_mm","exp_group_mm2","exp_group_mm3","exp_group_tert")]
surv_data_lin[,"exp_group_tert"]=as.vector(surv_data_lin[,"exp_group_tert"])
surv_data_lin[which(surv_data_lin[,"exp_group_tert"]=="low"),"exp_group_tert"]=1
surv_data_lin[which(surv_data_lin[,"exp_group_tert"]=="int"),"exp_group_tert"]=2
surv_data_lin[which(surv_data_lin[,"exp_group_tert"]=="high"),"exp_group_tert"]=3
surv_data_lin[,"exp_group_tert"]=as.numeric(surv_data_lin[,"exp_group_tert"])
survpvalue_linear=summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~exp_group_tert, data=surv_data_lin))$sctest[3]
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3)

#plot K-M estimate
pdf(file=survival_outfile4)
colors = rainbow(3)
title=paste("Survival by",gene,"expression (tertiles)")
plot(krfit.by_exp_group_tert, col = colors, xlab = "Time (years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)
groups=sort(unique(surv_data[,"exp_group_tert"])) #returns unique factor levels sorted alphabetically
names(colors)=groups
groups_custom=c("low","int","high")
colors_custom=colors[groups_custom]
group_sizes_custom=table(surv_data[,"exp_group_tert"])[groups_custom]
groups_custom=c("Low","Intermediate","High") #Reset names for consistency with manuscript
legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue_linear, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors_custom,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()


#Other tests to try
#coxph(Surv(general_t_rfs, general_e_rfs)~RF_predictions_responses, data=tamoxsens_data, method="breslow") #linear trend test?
#survdiff(tamoxsens_data.surv ~ RF_predictions_responses, data = tamoxsens_data, rho=1)
