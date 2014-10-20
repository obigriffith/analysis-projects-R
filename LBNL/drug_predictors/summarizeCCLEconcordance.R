#This script will summarize the concordance between drug response predictions from Gray lab predictors for an independent set of cell lines versus reported sensitivity from CCLE

#Load drug response metrics reported by CCLE
CCLE_drugdata_file="/Users/ogriffit/Dropbox/drug_predictors/CCLE/CCLE_NP24.2009_Drug_data_2012.02.20.csv"
CCLE_drugdata=read.csv(file=CCLE_drugdata_file)

#Set working dir for results
setwd("/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/Results/")
setwd("/Users/ogriffit/Dropbox/drug_predictors/CCLE/RtoolboxResults/Breast/Results_Fixed/")

#Load drug response predictions from Gray lab Rtoolbox
Rtoolbox_results_file="CCLE_results_ExpCNV.txt"
CCLE_predictions=read.table(file=Rtoolbox_results_file, header=TRUE, sep="\t")

#Drugs to analyze for concordance
#Graydrugs=c("17-AAG","AZD6244","Erlotinib","GSK_Tykerb","Nutlin 3a","Paclitaxel","PF-2341066","Sorafenib","TPT(FD)")
#CCLEdrugs=c("17-AAG","AZD6244","Erlotinib","Lapatinib","Nutlin-3","Paclitaxel","PF2341066","Sorafenib","Topotecan")

Graydrugs=c("Erlotinib","GSK_Tykerb","Paclitaxel","Sorafenib")
CCLEdrugs=c("Erlotinib","Lapatinib","Paclitaxel","Sorafenib")

pdf(file="CCLE_Concordance.pdf")
layout_setup = layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE), respect=TRUE)
#layout.show(layout_setup)

for (i in 1:length(Graydrugs)){
	Graydrug=Graydrugs[i]
	CCLEdrug=CCLEdrugs[i]
	#Limit to Rtoolbox (RTB) predictions and CCLE drug data for drug of interest
	RTB_pred_data=CCLE_predictions[which(CCLE_predictions[,"Drug.Compound"]==Graydrug),]
	CCLE_data=CCLE_drugdata[which(CCLE_drugdata[,"Compound"]==CCLEdrug),]
	rownames(CCLE_data)=CCLE_data[,"CCLE.Cell.Line.Name"]
	rownames(RTB_pred_data)=RTB_pred_data[,"Sample"]
	
	#Create waterfall plot for predicted sensitivity probability values for all cell lines tested in RTB (had CCLE molecular data but not necessarily CCLE drug data)
	cell_lines_all=as.vector(RTB_pred_data[,"Sample"])
	cell_line_names_all=sub("_BREAST","",cell_lines_all)
	RTB_preds_all=RTB_pred_data[cell_lines_all,"Probability"]

	cell_lines_all_sorted=cell_lines_all[order(RTB_preds_all, decreasing=TRUE)] #sort automatically excludes those with NA drug data
	cell_line_all_sorted_names=sub("_BREAST","",cell_lines_all_sorted)
	RTB_preds_all_sorted=sort(RTB_preds_all, decreasing=TRUE)
	names(RTB_preds_all_sorted)=cell_line_all_sorted_names
	ymin=floor(min(RTB_preds_all_sorted))
	ymax=ceiling(max(RTB_preds_all_sorted))
	barplot_values=barplot(RTB_preds_all_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.4, ylab="Sensitivity Probability", main=Graydrug)
	barplot(RTB_preds_all_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.4, ylab="Sensitivity Probability", main=Graydrug)
	box()	
	
	#Determine subset which actually have drug data from CCLE to correlate with
	cell_lines=as.vector(RTB_pred_data[which(RTB_pred_data[,"Sample"]%in%CCLE_data[,"CCLE.Cell.Line.Name"]),"Sample"])
	cell_line_names=sub("_BREAST","",cell_lines)
	
	CCLE_EC50=CCLE_data[cell_lines,"EC50..uM."]
	CCLE_EC50_neglog=-log10(CCLE_EC50)
	CCLE_IC50=CCLE_data[cell_lines,"IC50..uM."]
	CCLE_IC50_neglog=-log10(CCLE_IC50)
	RTB_preds=RTB_pred_data[cell_lines,"Probability"]

	#Create corresponding waterfall plot for predicted sensitivity probability values
	cell_lines_sorted=cell_lines[order(RTB_preds, decreasing=TRUE)] #sort automatically excludes those with NA drug data
	cell_line_sorted_names=sub("_BREAST","",cell_lines_sorted)
	RTB_preds_sorted=sort(RTB_preds, decreasing=TRUE)
	names(RTB_preds_sorted)=cell_line_sorted_names
	ymin=floor(min(RTB_preds_sorted))
	ymax=ceiling(max(RTB_preds_sorted))
	barplot_values=barplot(RTB_preds_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="Sensitivity Probability", main=Graydrug)
	barplot(RTB_preds_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="Sensitivity Probability", main=Graydrug)
	box()

	#Create corresponding waterfall plot for IC50 values
	cell_lines_sorted=cell_lines[order(CCLE_IC50_neglog, decreasing=TRUE)] #sort automatically excludes those with NA drug data
	cell_line_sorted_names=sub("_BREAST","",cell_lines_sorted)
	CCLE_IC50_neglog_sorted=sort(CCLE_IC50_neglog, decreasing=TRUE)
	names(CCLE_IC50_neglog_sorted)=cell_line_sorted_names
	ymin=floor(min(CCLE_IC50_neglog_sorted))
	ymax=ceiling(max(CCLE_IC50_neglog_sorted))
	barplot_values=barplot(CCLE_IC50_neglog_sorted, plot=FALSE, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="-log10(IC50)", main=Graydrug)
	barplot(CCLE_IC50_neglog_sorted, ylim=c(ymin,ymax), xpd=FALSE, las=2, cex.names=0.7, ylab="-log10(IC50)", main=Graydrug)
	box()

	#Create correlation plot
	#print(cbind(CCLE_EC50,RTB_preds))
	#plot(x=CCLE_EC50,y=RTB_preds, main=Graydrug)
	plot(x=CCLE_IC50_neglog,y=RTB_preds, main=Graydrug, ylab="Sensitivity Probability", xlab="-log10(CCLE IC50)")
	r=cor(x=CCLE_IC50_neglog,y=RTB_preds,method="spearman")
	legend("top", legend=paste("r=",format(r,digits=4),sep=""), bty="n")
}
dev.off()





