
#Extract data for specific gene lists from specific data types 
#Compare resistant vs sensitive for specific drug

#results file
results_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ClueGO/NegRegKinase_pathway_gene_analysis.txt"

#Set path to data files
exon_array_gene_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/ExonArray/filtered_v2/breastExon_genelevel_stringent.csv"
meth_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/methylation/filtered_v2/Methylation_stringent.csv"
meth_annotation_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/methylation/filtered_v2/Methylation_annotation_stringent.csv"
RNAseq_gene_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RNAseq/filtered_v3/breastRNAseq_genelevel_stringent.txt"
SNP6_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/CNV/filtered_v2/SNP6_genelevel_stringent_std0.7.csv"
U133A_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/U133A/filtered_v2/Neve_AffyRMA_genelevel_maxvar_stringent.csv"

#import molecular data
CellLineEAData=read.csv(exon_array_gene_file,sep=",",header=TRUE, row.names=1)
CellLineU133AData=read.csv(U133A_file,sep=",",header=TRUE, row.names=1)
CellLineRSData=read.table(RNAseq_gene_file,sep="\t",header=TRUE, row.names=3)
CellLineRSData=CellLineRSData[,4:length(colnames(CellLineRSData))]
CellLineCNVData=read.csv(SNP6_file,sep=",",header=TRUE,row.names=5)
CellLineCNVData=CellLineCNVData[,5:length(colnames(CellLineCNVData))]
CellLineMethData=read.csv(meth_file,sep=",",header=TRUE, row.names=1)
GeneProbeMapping.meth=read.csv(meth_annotation_file,sep=",",header=TRUE)
CellLineMethData.aggr=aggregate.data.frame(CellLineMethData, by=list(GeneProbeMapping.meth[,"Symbol"]), mean)
rownames(CellLineMethData.aggr)=CellLineMethData.aggr[,"Group.1"]
CellLineMethData=CellLineMethData.aggr[,-1]

#Fix column names
colnames(CellLineEAData)[which(colnames(CellLineEAData)=="X184A1")]="184A1"
colnames(CellLineEAData)[which(colnames(CellLineEAData)=="X184B5")]="184B5"
colnames(CellLineEAData)[which(colnames(CellLineEAData)=="X600MPE")]="600MPE"
colnames(CellLineU133AData)[which(colnames(CellLineU133AData)=="X600MPE")]="600MPE"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X184A1")]="184A1"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X184B5")]="184B5"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X600MPE")]="600MPE"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X21MT1")]="21MT1"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X21MT2")]="21MT2"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X21NT")]="21NT"
colnames(CellLineRSData)[which(colnames(CellLineRSData)=="X21PT")]="21PT"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X184A1")]="184A1"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X184B5")]="184B5"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X600MPE")]="600MPE"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X21MT1")]="21MT1"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X21MT2")]="21MT2"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X21NT")]="21NT"
colnames(CellLineCNVData)[which(colnames(CellLineCNVData)=="X21PT")]="21PT"
colnames(CellLineMethData)[which(colnames(CellLineMethData)=="X600MPE")]="600MPE"

#Import drug cutoff data
#Filter down to core set of cell lines - must have drug data and at least one other molecular profiling data type
core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")

#drug response data
drugdatafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/LBNL_Gray_BCCL_alldrugs_26Aug11.csv"
raw_drugdata_import=read.csv(drugdatafile)

drug_data=raw_drugdata_import[,2:length(colnames(raw_drugdata_import))]
rownames(drug_data)=as.vector(raw_drugdata_import[,1])

#Retrieve data for only libraries in core cell line set
drug_data_filt=drug_data[core_cell_lines,]

#Transform data to -log10 values
drug_data_filt_trans=-log10(drug_data_filt)
drugs=colnames(drug_data_filt_trans)

#Retrieve mean GI50 values
GI50meanfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/drugdata/GI50meanThresholds_v2.csv"
GI50mean_import=read.csv(GI50meanfile)
GI50mean_data=GI50mean_import[,2]
names(GI50mean_data)=GI50mean_import[,1]
drug_names=GI50mean_import[,1]
cbind(drugs,names(GI50mean_data)) #Check to make sure drugs are in same order in drug data and mean GI50 files


#Pathway="negative regulation of protein kinase"
#See ~/Dropbox/drug_predictors/ClueGOPathwaysDrugCounts_noRPPA.xlsx for pathway and drugs associated with
#See ~/Dropbox/drug_predictors/PathwayAnalysis_V2.docx for predictor model/signature behind pathway analysis for each drug
#See ~/Dropbox/Manuscript_drug_predictors/NatureSupplementary/SupplTable10.xlsx for actual genes from this pathway in each signature

##AKT1-2 inhibitor (Exon array (gene level), RF)##
path_genes=c("CAV1","DEPDC6","DUSP7","ERRFI1","GADD45A","PDCD4","PRKCA","SFRP1","SPRY1","SPRY2","UCHL1")
drug="AKT1.2.inhibitor"
drug_name="AKT1-2 inhibitor"
data_type="Exon array (gene level)"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]
resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

path_gene_data=CellLineEAData[path_genes,]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results1=data.frame(cbind(path_genes=path_genes, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes)){
 path_gene=path_genes[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results1[i,"drug"]=drug_name
 results1[i,"datatype"]=data_type
 results1[i,"mean_sens"]=mean_sensitive
 results1[i,"mean_res"]=mean_resistant
 results1[i,"FC"]=FC
 results1[i,"pvalue"]=wilcox.pval
}


##Bosutinib (U133A + SNP6, RF)##
path_genes=c("DUSP6","GADD45A","PDCD4","SFRP1")
drug="SKI.606.Bosutinib."
drug_name="SKI-606(Bosutinib)"
data_type="U133A"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]

##U133A
resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineU133AData[path_genes,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results2=data.frame(cbind(path_genes=path_genes, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes)){
 path_gene=path_genes[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results2[i,"drug"]=drug_name
 results2[i,"datatype"]=data_type
 results2[i,"mean_sens"]=mean_sensitive
 results2[i,"mean_res"]=mean_resistant
 results2[i,"FC"]=FC
 results2[i,"pvalue"]=wilcox.pval
}

##SNP6 - No genes found in SNP6 data
CellLineCNVData[path_genes,]


##TCS PIM-11 (U133A + methylation, RF)
path_genes=c("CAV1","CDKN2A","DAB2IP","DBNDD2","DRD1","DUSP22","DUSP8","ERRFI1","GADD45G","HEXIM1","LAX1","MAPK8IP1","NF2","NR2F2","PAK3","PDCD4","PDPK1","PKIG","PRKAG2","PSEN1","RB1","RGS3","RGS4","SFN","SFRP1","SPRY1","TAF7","TARBP2","TP73","TRIB2","TRIB3","UCHL1","ZGPAT")
drug="TCS.PIM.11"
drug_name="TCS PIM-11"
data_type="U133A"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]

##U133A
path_genes_U133A=path_genes[which(path_genes %in% rownames(CellLineU133AData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineU133AData[path_genes_U133A,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results3=data.frame(cbind(path_genes=path_genes_U133A, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_U133A)){
 path_gene=path_genes_U133A[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results3[i,"drug"]=drug_name
 results3[i,"datatype"]=data_type
 results3[i,"mean_sens"]=mean_sensitive
 results3[i,"mean_res"]=mean_resistant
 results3[i,"FC"]=FC
 results3[i,"pvalue"]=wilcox.pval
}

##Meth
data_type="Meth"
path_genes_Meth=path_genes[which(path_genes %in% rownames(CellLineMethData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineMethData[path_genes_Meth,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results4=data.frame(cbind(path_genes=path_genes_Meth, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_Meth)){
 path_gene=path_genes_Meth[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results4[i,"drug"]=drug_name
 results4[i,"datatype"]=data_type
 results4[i,"mean_sens"]=mean_sensitive
 results4[i,"mean_res"]=mean_resistant
 results4[i,"FC"]=FC
 results4[i,"pvalue"]=wilcox.pval
}


##GSK2119563 (Exon array (gene level), RF)##
path_genes=c("CAV1","PDCD4","SFRP1","SPRY1","SPRY2","SPRY4","TRIB2","UCHL1")
drug="GSK2119563A"
drug_name="GSK2119563A"
data_type="Exon array (gene level)"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]
resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

path_gene_data=CellLineEAData[path_genes,]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results5=data.frame(cbind(path_genes=path_genes, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes)){
 path_gene=path_genes[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results5[i,"drug"]=drug_name
 results5[i,"datatype"]=data_type
 results5[i,"mean_sens"]=mean_sensitive
 results5[i,"mean_res"]=mean_resistant
 results5[i,"FC"]=FC
 results5[i,"pvalue"]=wilcox.pval
}


##GSK2126458 (Exon array (gene level), RF)##
path_genes=c("BMP4","CAV1","DEPDC6","DUSP7","GADD45A","PDCD4","PRKCA","SFRP1","SPRY1","SPRY2","TRIB2","UCHL1")
drug="GSK2126458A"
drug_name="GSK2126458A"
data_type="Exon array (gene level)"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]
resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

path_gene_data=CellLineEAData[path_genes,]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results6=data.frame(cbind(path_genes=path_genes, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes)){
 path_gene=path_genes[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results6[i,"drug"]=drug_name
 results6[i,"datatype"]=data_type
 results6[i,"mean_sens"]=mean_sensitive
 results6[i,"mean_res"]=mean_resistant
 results6[i,"FC"]=FC
 results6[i,"pvalue"]=wilcox.pval
}


##PF-4691502 (Exon array (gene level) + meth, LS-SVM)##
path_genes=c("APOE","CAV1","DEPDC6","DRD1","DUSP7","GADD45A","PDCD4","SPRY1")
drug="PF.4691502"
drug_name="PF-4691502"
data_type="Exon array (gene level)"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]

##exon-array
path_genes_EA=path_genes[which(path_genes %in% rownames(CellLineEAData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineEAData[path_genes_EA,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results7=data.frame(cbind(path_genes=path_genes_EA, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_EA)){
 path_gene=path_genes_EA[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results7[i,"drug"]=drug_name
 results7[i,"datatype"]=data_type
 results7[i,"mean_sens"]=mean_sensitive
 results7[i,"mean_res"]=mean_resistant
 results7[i,"FC"]=FC
 results7[i,"pvalue"]=wilcox.pval
}

##Meth
data_type="Meth"
path_genes_Meth=path_genes[which(path_genes %in% rownames(CellLineMethData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineMethData[path_genes_Meth,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results8=data.frame(cbind(path_genes=path_genes_Meth, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_Meth)){
 path_gene=path_genes_Meth[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results8[i,"drug"]=drug_name
 results8[i,"datatype"]=data_type
 results8[i,"mean_sens"]=mean_sensitive
 results8[i,"mean_res"]=mean_resistant
 results8[i,"FC"]=FC
 results8[i,"pvalue"]=wilcox.pval
}


##TGX-221 (RNA-seq (gene level) + meth + SNP6, RF)##
path_genes=c("APOE","DBNDD2","IL1B","RGS4")
drug="GSK1059868A"
drug_name="GSK1059868A"
data_type="RNA-seq (gene level)"

#Divide cell lines into responders and non-responders
drug_data_interest=as.numeric(drug_data_filt_trans[,drug])
names(drug_data_interest)=rownames(drug_data_filt_trans)
cell_lines_sorted=names(sort(drug_data_interest, decreasing=TRUE)) #sort automatically excludes those with NA drug data
drug_data_interest_sorted=drug_data_interest[cell_lines_sorted]
mean_cutoff=GI50mean_data[drug_name]

##RNA-seq
path_genes_RS=path_genes[which(path_genes %in% rownames(CellLineRSData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineRSData[path_genes_RS,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results9=data.frame(cbind(path_genes=path_genes_RS, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_RS)){
 path_gene=path_genes_RS[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results9[i,"drug"]=drug_name
 results9[i,"datatype"]=data_type
 results9[i,"mean_sens"]=mean_sensitive
 results9[i,"mean_res"]=mean_resistant
 results9[i,"FC"]=FC
 results9[i,"pvalue"]=wilcox.pval
}

##Meth
data_type="Meth"
path_genes_Meth=path_genes[which(path_genes %in% rownames(CellLineMethData))]

resistants=names(which(drug_data_interest<=mean_cutoff))
sensitives=names(which(drug_data_interest>mean_cutoff))

#Extract data for analysis and determine cell lines in that datatype
path_gene_data=CellLineMethData[path_genes_Meth,]
sensitives=sensitives[which(sensitives %in% colnames(path_gene_data))]
resistants=resistants[which(resistants %in% colnames(path_gene_data))]

#Test each gene for difference between sensitives and resistants
#Create data frame to store result
results10=data.frame(cbind(path_genes=path_genes_Meth, drug=NA, datatype=NA, mean_sens=NA, mean_res=NA, FC=NA, pvalue=NA), stringsAsFactors=FALSE)

for (i in 1:length(path_genes_Meth)){
 path_gene=path_genes_Meth[i]
 mean_resistant=mean(as.numeric(path_gene_data[path_gene,resistants]))
 mean_sensitive=mean(as.numeric(path_gene_data[path_gene,sensitives]))
 FC=mean_resistant/mean_sensitive
 if(FC<1){FC = -1/FC}
 wilcox.pval=wilcox.test(x=as.numeric(path_gene_data[path_gene,resistants]), y=as.numeric(path_gene_data[path_gene,sensitives]))$p.value
 results10[i,"drug"]=drug_name
 results10[i,"datatype"]=data_type
 results10[i,"mean_sens"]=mean_sensitive
 results10[i,"mean_res"]=mean_resistant
 results10[i,"FC"]=FC
 results10[i,"pvalue"]=wilcox.pval
}

final_results=rbind(results1,results2,results3,results4,results5,results6,results7,results8,results9,results10)
write.table(final_results, file=results_file, row.names=FALSE, sep="\t")
