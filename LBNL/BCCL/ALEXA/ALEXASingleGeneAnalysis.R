#Cell line info
celllinefile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/cell_line_info.2.txt"

#Gene data
genedatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/MAP3K14/ENSG00000006062.txt"

#All gene expression data
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"

#Set output dir for feature type
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/MAP3K14/"
setwd(outdir)

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_cell_data=read.table(celllinefile, header = TRUE, na.strings = "NA", sep="\t")
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_gene_data=read.table(genedatafile, header = TRUE, na.strings = "NA", sep="\t")

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:length(colnames(raw_data_import))]
raw_exp_status=raw_exp_status_import[,5:length(colnames(raw_data_import))]

#Make sure that cell line info and raw data are ordered the same!!!
libs=colnames(raw_data)
lib_names=as.vector(raw_cell_data[,"Sample.Name"])
cbind(lib_names,libs)


#All libraries
data=raw_data
feat_data=raw_feat_data
exp_status=raw_exp_status
cell_data=raw_cell_data

#Exclude low quality (and "Normal" and "Other" cell lines)
#high_qual=which(as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int")
high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"Sample.Type"]!="Other") & as.vector(raw_cell_data[,"subtype"]!="Normal"))
#high_qual=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"Sample.Type"]!="Other") & as.vector(raw_cell_data[,"subtype"]!="Normal") & as.vector(raw_cell_data[,"subtype"]!="Basal_NM"))

#Exclude low quality and "Unknown/Normal subtype" libs
high_qual_known=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"subtype"]!="Unknown") & as.vector(raw_cell_data[,"subtype"]!="Normal"))

#Exclude low quality and "Unknown/Normal/Normal-like subtype" libs
#high_qual_known=which((as.vector(raw_cell_data[,"Quality"])=="High" | as.vector(raw_cell_data[,"Quality"])=="Int") & as.vector(raw_cell_data[,"subtype"]!="Unknown") & as.vector(raw_cell_data[,"subtype"]!="Normal") & as.vector(raw_cell_data[,"subtype"]!="Basal_NM"))

#Apply library filter to datasets
data_HQK=raw_data[,high_qual_known]
feat_data_HQK=raw_feat_data
exp_status_HQK=raw_exp_status[,high_qual_known]
cell_data_HQK=raw_cell_data[high_qual_known,]
libs=colnames(data_HQK)

#Retrieve cell line details
lib_names=as.vector(cell_data_HQK[,"Sample.Name"])
subtypes=as.vector(cell_data_HQK[,"subtype"])
ERBB2=as.vector(cell_data_HQK[,"ERBB2"])
Qualities=as.vector(cell_data_HQK[,"Quality"])
cbind(lib_names,libs)

#Set up colors for sidebars, etc
colors=rainbow(6)
#Subtype
subtype_colors=subtypes
subtype_colors[subtype_colors=="Basal"]=colors[1] #red
subtype_colors[subtype_colors=="Luminal"]=colors[2] #yellow
subtype_colors[subtype_colors=="ClaudinLow"]=colors[3] #green
subtype_colors[subtype_colors=="Basal_NM"]=colors[4] #blue
subtype_colors[subtype_colors=="Normal"]=colors[5]
subtype_colors[subtype_colors=="Unknown"]=colors[6]

#HER2 status
her2_colors=ERBB2
her2_colors[her2_colors=="HighAmp"]=colors[1]
her2_colors[her2_colors=="OverExp"]=colors[2]
her2_colors[her2_colors=="LowAmp"]=colors[3]
her2_colors[her2_colors=="NoAmp"]=colors[5]
her2_colors[her2_colors=="ND"]=colors[6]

#Quality status
quality_colors=Qualities
quality_colors[quality_colors=="High"]=colors[1]
quality_colors[quality_colors=="Int"]=colors[3]
quality_colors[quality_colors=="Low"]=colors[5]

#color mapping for later plots
colormap=c(colors[1],colors[2],colors[3],colors[4])
names(colormap)=c("Basal","Luminal","ClaudinLow","Basal_NM")
subtype_names=c("Basal","Luminal","Claudin-Low","Normal-Like")
names(subtype_names)=c("Basal","Luminal","ClaudinLow","Basal_NM")

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
pe_data=apply(exp_status_HQK, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data_filt=data_HQK[passed_pe,]
feat_data_filt=feat_data_HQK[passed_pe,]

cov_data=apply(data_filt, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data_filt=data_filt[passed_cov,]
feat_data_filt=feat_data_filt[passed_cov,]
#rownames(data_filt)=feat_data_filt[,"Seq_Name"]
rownames(data_filt)=feat_data_filt[,"EnsEMBL_Gene_ID"] #Use for gene-level data because SeqName not unique

#Convert values to log2, unless they are 0 in which case set them to 0
z = log2(data_filt+1)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Heatmap
pdf(file=heatmap_outfile)
x=t(as.matrix(z))
#x=t(x)
main_title="RNAseq gene-level expression"
#main_title="RNAseq transcript-level expression"
rlab=cbind(subtype_colors,her2_colors,quality_colors)
colnames(rlab) = c("Subtype","Her2","Quality")
heatmap.plus(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", main=main_title, RowSideColors=rlab, labRow=lib_names, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
dev.off()

pdf(file=heatmap_outfile2)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, RowSideColors=subtype_colors, labRow=lib_names, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
dev.off()

#attempt to classify by subtype
#First exclude some cell lines that do not cluster where expected above. We want to train classifier on only those lines we are confident about
data_train=z[,!colnames(z)%in%c("BT20","HCC1954","MDAMB157","SUM159PT","MDAMB175VII")]
subtypes_train=subtypes[!colnames(z)%in%c("BT20","HCC1954","MDAMB157","SUM159PT","MDAMB175VII")]
lib_names_train=colnames(data_train)

subtypes_unique=unique(subtypes)

#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(subtypes_train))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

#rf_model=randomForest(x=t(data_train), y=as.factor(subtypes_train), importance = TRUE, ntree = 50001, proximity=TRUE)
rf_model=randomForest(x=t(data_train), y=as.factor(subtypes_train), importance = TRUE, ntree = 50001, proximity=TRUE, sampsize=sampsizes)

#Save model for later use:
save(rf_model, file=rf_model_file)

#Get importance measures
rf_importances=importance(rf_model, scale=FALSE)

#Determine performance statistics
confusion=rf_model$confusion
overall_error=rf_model$err.rate[length(rf_model$err.rate[,1]),1]*100
overall_accuracy=100-overall_error

#Prepare stats for output to file
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")

#Prepare confusion table for writing to file
confusion_out=confusion[1:num_classes,1:num_classes]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,num_classes+2],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(acc_out,err_out), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_model, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
target_labels=as.vector(subtypes_train)
target_labels[target_labels=="Basal_NM"]="N"
target_labels[target_labels=="Basal"]="B"
target_labels[target_labels=="Luminal"]="L"
target_labels[target_labels=="ClaudinLow"]="C"
pdf(file=MDS_pdffile)
MDSplot(rf_model, as.factor(subtypes_train), k=2, xlab="", ylab="", pch=target_labels, palette=colormap[sort(unique(subtypes_train))], main="MDS plot")
dev.off()

#plot distributions of scores for each predicted subtype
pdf(file=vote_dist_pdffile)
subtype_color=colormap[sort(unique(subtypes_train))]
par(mfrow=c(num_classes,num_classes))
for (i in 1:num_classes){
  true_subtype=rf_model$classes[i]
  for (j in 1:num_classes){
    pred_subtype=rf_model$classes[j]
    if(j==1){ylabel=subtype_names[true_subtype]}else{ylabel=""}
    if(i==num_classes){xlabel=paste(subtype_names[pred_subtype],"Score")}else{xlabel=""}
    hist(rf_model$votes[rf_model$y==true_subtype,pred_subtype],xlim=c(0,1),main="", xlab=xlabel, ylab=ylabel, col=subtype_color[j])
  }
}
dev.off()

#Create heatmap of just the training samples for top XXX genes
#Subtype
subtype_colors_train=subtypes_train
subtype_colors_train[subtype_colors_train=="Basal"]=colors[1]
subtype_colors_train[subtype_colors_train=="Luminal"]=colors[2]
subtype_colors_train[subtype_colors_train=="ClaudinLow"]=colors[3]
subtype_colors_train[subtype_colors_train=="Basal_NM"]=colors[4]
subtype_colors_train[subtype_colors_train=="Normal"]=colors[5]
subtype_colors_train[subtype_colors_train=="Unknown"]=colors[6]

#Limit to top XXX genes - 300 determined to be approx optimal with experimentation for producing four clean clusters with min number of genes
topgenes=sort(rf_importances[,num_classes+2], decreasing=TRUE)[1:300]
data_train_topgenes=data_train[names(topgenes),]

pdf(file=heatmap_outfile3)
x=t(as.matrix(data_train_topgenes))
main_title="RNAseq gene-level expression, training, top 300 genes"
par(cex.main=0.9)
heatmap.2(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, RowSideColors=subtype_colors_train, labRow=lib_names_train, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
dev.off()

#Apply RF subtype predictor to all lines including problem cell lines and Unknowns, exclude low quality 
#Note: predicted classes and probability should be used from above for training data and then combined with results for test data here

#Apply library filter to datasets
data_HQ=raw_data[,high_qual]
feat_data_HQ=raw_feat_data
exp_status_HQ=raw_exp_status[,high_qual]
cell_data_HQ=raw_cell_data[high_qual,]
libs_HQ=colnames(data_HQ)
rownames(data_HQ)=feat_data_HQ[,"EnsEMBL_Gene_ID"] #Use for gene-level data because SeqName not unique

#Retrieve cell line details
lib_names_HQ=as.vector(cell_data_HQ[,"Sample.Name"])
subtypes_HQ=as.vector(cell_data_HQ[,"subtype"])
ERBB2_HQ=as.vector(cell_data_HQ[,"ERBB2"])
Qualities_HQ=as.vector(cell_data_HQ[,"Quality"])
cbind(lib_names_HQ,libs_HQ)

#Run data through forest
#Convert values to log2, unless they are 0 in which case set them to 0
predictor_data = t(log2(data_HQ+1))
RF_predictions_response=predict(rf_model, predictor_data, type="response")
RF_predictions_prob=predict(rf_model, predictor_data, type="prob")
RF_predictions_vote=predict(rf_model, predictor_data, type="vote", norm.votes=FALSE)

#Combine predictions from training step and test step
all_lines=colnames(data_HQ)
train_lines=colnames(data_train)
test_lines=all_lines[!all_lines%in%train_lines]
test_lines_ind=which(!all_lines%in%train_lines)

#Combine train/test predictions and then put back into original order
RF_combined_prob=rbind(rf_model$votes,RF_predictions_prob[test_lines,])
RF_combined_response=c(as.vector(rf_model$y),as.vector(RF_predictions_response[test_lines_ind]))
RF_combined=cbind(RF_combined_prob,RF_combined_response)[all_lines,] #Combine into one dataframe to allow sorting of probs and responses correctly
RF_combined_prob=apply(RF_combined[,1:num_classes],2,as.numeric) #Now, that sorted, pull out separately again and reformat as before
rownames(RF_combined_prob)=rownames(RF_combined)
RF_combined_response=as.factor(RF_combined[,num_classes+1]) #Now, that sorted, pull out separately again and reformat as before

#Plot probabilities for each library
pdf(file=prob_pdffile, height=7, width=10)
par(oma=c(2,1,1,1), lab=c(5,10,7))#oma = bottom, left, top, right
order=order(x= -xtfrm(RF_combined_response), y=apply(RF_combined_prob,1,max), decreasing=TRUE)
barplot(t(RF_combined_prob[order,]), col=colormap[sort(unique(subtypes_train))], las=3, cex.names=0.8, ylab="Probability")
grid(nx = NA, ny = NULL, col = "black")
legend(x=50, y=1, legend=subtype_names, fill=colormap, bg="white", border="white")
dev.off()

#Add concept of uncertainty. Don't allow assignment to subtype unless probability for that subtype is at least > some threshold
#Re assign such cases to unknown for now?
#Two methods:
#1. Determine difference between highest and next-highest probability
best_prob_diff_fun=function(x){
  best_prob=max(x)
  next_best_prob=sort(x)[length(x)-1]
  best_prob_diff=best_prob - next_best_prob
  weighted_prob=best_prob * (1-next_best_prob)
  return(best_prob_diff)
}
#2. Determine a weighted probability where best probability is penalized by next-highest probability
weighted_prob_fun=function(x){
  best_prob=max(x)
  next_best_prob=sort(x)[length(x)-1]
  weighted_prob=best_prob * (1-next_best_prob)
  return(weighted_prob)
}
best_prob_diff=apply(RF_combined_prob, 1, best_prob_diff_fun)
weighted_prob=apply(RF_combined_prob, 1, weighted_prob_fun)
uncertain=names(which(weighted_prob<0.28)) #Threshold of 0.28 chosen by manual inspection

#Add best_prob_diff and weighted_prob to dataframe and reassign uncertain classes to "Unknown"
RF_combined=cbind(RF_combined,best_prob_diff,weighted_prob)
RF_combined[uncertain,"RF_combined_response"]="Unknown"

#Create new heatmap showing all cell lines (except low qual) with prior and predicted subtypes
#Limit to top X genes
topgenes=sort(rf_importances[,num_classes+2], decreasing=TRUE)[1:500]

#Subtype colors - prior
subtype_colors_prior=subtypes_HQ
subtype_colors_prior[subtype_colors_prior=="Basal"]=colors[1]
subtype_colors_prior[subtype_colors_prior=="Luminal"]=colors[2]
subtype_colors_prior[subtype_colors_prior=="ClaudinLow"]=colors[3]
subtype_colors_prior[subtype_colors_prior=="Basal_NM"]=colors[4]
subtype_colors_prior[subtype_colors_prior=="Normal"]=colors[5]
subtype_colors_prior[subtype_colors_prior=="Unknown"]=colors[6]

#Subtype colors - predicted
subtype_colors_pred=as.vector(RF_combined[,"RF_combined_response"])
subtype_colors_pred[subtype_colors_pred=="Basal"]=colors[1]
subtype_colors_pred[subtype_colors_pred=="Luminal"]=colors[2]
subtype_colors_pred[subtype_colors_pred=="ClaudinLow"]=colors[3]
subtype_colors_pred[subtype_colors_pred=="Basal_NM"]=colors[4]
subtype_colors_pred[subtype_colors_pred=="Normal"]=colors[5]
subtype_colors_pred[subtype_colors_pred=="Unknown"]=colors[6]

pdf(file=heatmap_outfile4)
x=as.matrix(predictor_data[,names(topgenes)])
rlab=cbind(subtype_colors_prior,subtype_colors_pred)
colnames(rlab) = c("Prior","Predicted")
main_title="RNAseq gene-level expression, predicted subtypes"
par(cex.main=1)
heatmap.plus(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", main=main_title, RowSideColors=rlab, labRow=lib_names_HQ, labCol=FALSE, cexRow=0.8, cexCol=0.60, col=rev(heat.colors(75)))
legend("topleft", legend=c("Basal","Luminal","Claudin-low","Normal-like","Unknown"), fill=c(colors[1],colors[2],colors[3],colors[4],colors[6]), bty="n", cex=0.6, border="white")
dev.off()

#Print results to file
RF_combined=cbind(rownames(RF_combined),RF_combined)
write.table(RF_combined,file=subtype_pred_outfile, sep="\t", quote=FALSE, row.names=FALSE)

