library(randomForest)
library(gregmisc)

#Set working directory and filenames for Input/output
setwd("C:/Documents and Settings/obig/My Documents/Projects/Thyroid_markers/BenignMalignantCombined/benign_malignant_paper/updated_analysis/data_files")
datafile="BenignMalignant_58markers_23JAN08_ungrouped_and_grouped.txt"
#outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput_2marker_combos.txt"
#outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput_3marker_combos.txt"
#outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput_4marker_combos.txt"
#outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput_2marker_combos_top7.txt"
outfile="BenignMalignant_58markers_23JAN08_ungrouped_NoM_RFoutput_3marker_combos_top7.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Set the number of genes to try for each combination (e.g., 2-gene panel, 3-gene panel, etc)
#n_gene=2
n_gene=3
#n_gene=4

#Exclude any cases you do not want used in the classifier (e.g. in this case we do not want any M or HCC tumors where path=3 or 4)
data=data[data[,13]<4,] #Keep only rows where pathology is 0, 1, 2 or 3

#Get potential predictor variables from datafile
#predictor_data=data[,36:92]
#predictor_data=data[,c("Galectin","CK19","VEGF","AuroraA","P16","AR","HBME","BCL2","CYCLIND1","Caveolin1","CYCLINE","ECAD","CR3","Clusterin","IGFBP5","P21","IGFBP2","BetaCatenin","HER4","TG","KI67","Caveolin","AuroraC","S100","MRAS","CKIT","HER3","RET","AMFR","MLH1","AAT","TTF1","PGI","HSP27","Syntrophin")]
predictor_data=data[,c("AR","AuroraA","CK19","Galectin","HBME","P16","VEGF")]

#predictor_names=c("AAT","AMF-R","AR","Aurora-A","Aurora-C","Bcl-2","CTNNB1","CDX2","CK19","c-kit","COX2","CR3","Cyclin-D1","Cyclin-E","Caveolin","CAV-1","Clusterin","E-CAD","ER","Galectin-3","HBME-1","EGFR","HER2","HER3","HER4","HSP-27","IGFBP2","IGFBP5","INH","KI67","MDM2","MLH1","MRAS","O13","P16","P21","P27","P53","P57","P63","P75-NTR","PGI","PMS2","PR","PSA","RET","S100","Syntrophin","TDT","TG","TOPO-II","TS106","TSH","TTF-1","VEGF","WT1","P504S")
#predictor_names=c("Galectin-3","CK19","VEGF","Aurora-A","P16","AR","HBME-1","Bcl-2","Cyclin-D1","CAV-1","Cyclin-E","E-CAD","CR3","Clusterin","IGFBP5","P21","IGFBP2","CTNNB1","HER4","TG","KI67","Caveolin","Aurora-C","S100","MRAS","c-kit","HER3","RET","AMF-R","MLH1","AAT","TTF-1","PGI","HSP-27","Syntrophin")
predictor_names=c("AR","Aurora-A","CK19","Galectin-3","HBME-1","P16","VEGF")

colnames(predictor_data)=predictor_names

#Get target variable and specify as factor/categorical
target=(data[,14])
target[target==0]="Benign"
target[target==1]="Malignant"
target=as.factor(target)

#If there are predictor variables that are constant/invariant, consider removing them
#Make sure predictor data is correctly recognized as categorical and determine which variables are invariant (have only one value)
variant_list=vector(mode = "logical", length = 0)
for(i in 1:length(predictor_data)){
 predictor_data[,i]=as.factor(predictor_data[,i])
 if(nlevels(predictor_data[,i])<2){
 variant_list[i]=FALSE
 }
 if(nlevels(predictor_data[,i])>=2){
 variant_list[i]=TRUE
 }
}
predictor_data_variant=predictor_data[,variant_list]

#Get list of predictors still remaining (i.e., not invariant)
predictor_list=colnames(predictor_data_variant)

#Determine all n-gene combinations
combos=combinations(length(predictor_list), n_gene, predictor_list)

#Create list of gene-combination names (e.g., "AR/CK19/P16")
combo_names=vector(length=length(combos[,1]))
for (i in 1:length(combo_names)){
 #combo_names[i]=paste(combos[i,1],combos[i,2],combos[i,3],sep="/")
 combo_names[i]=paste(combos[i,1:length(combos[i,])], collapse="/")
}

#Create array to store results for each RF model
RF_results = array(0, dimnames = list(combo_names, c("Sensitivity", "Specificity", "PPV", "NPV", "accuracy", "Mis_Ben", "Mis_Mal")), dim=c(length(combo_names),7))

for(i in 1:length(combos[,1])){
 combo_predictor_data=predictor_data_variant[,combos[i,]]
 #If there are missing values these will have to be dealt with somehow. Use rfImpute 
 #Warning: rfImpute adds the target/response variable as the first column, therefore remove this column with "[,-1]"

 NA_check=is.na(combo_predictor_data)
 if(length(NA_check[NA_check==TRUE])>0){
  RF_predictor_data=rfImpute(x=combo_predictor_data, y=target, iter=5)[,-1]
 }else{
  RF_predictor_data=combo_predictor_data #If no NAs, just use data as is
 }

 #Run RandomForests with imputed marker data as potential predictors and with path_gp1 as target
 rf_output=randomForest(x=RF_predictor_data, y=target, importance = TRUE, ntree = 1000, proximity=TRUE)

 #Determine performance statistics
 confusion=rf_output$confusion
 RF_results[i,"Sensitivity"]=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
 RF_results[i,"Specificity"]=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
 RF_results[i,"PPV"]=(confusion[2,2]/(confusion[2,2]+confusion[1,2]))*100
 RF_results[i,"NPV"]=(confusion[1,1]/(confusion[1,1]+confusion[2,1]))*100
 overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
 RF_results[i,"accuracy"]=100-overall_error
 RF_results[i,"Mis_Ben"]=confusion[1,2]
 RF_results[i,"Mis_Mal"]=confusion[2,1]
}

#Write complete results to file
write.table (RF_results, sep="\t", file=outfile)
