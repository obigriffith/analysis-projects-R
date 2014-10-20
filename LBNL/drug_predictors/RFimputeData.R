library(randomForest)

outdir="C:/Users/Obi/Documents/My Dropbox/drug_predictors/combined/"
setwd(outdir)
outfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/combined/BCCL_combined_data.4.imp.txt"

#Import cell line data
cell_line_datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/BCCL_data_list_v2.txt"
raw_cell_line_import=read.table(cell_line_datafile, header = TRUE, na.strings = "NA", sep="\t")
rownames(raw_cell_line_import)=raw_cell_line_import[,1]

#Import combined predictor data
datafile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/combined/BCCL_combined_data.4.txt"
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:2]
raw_data=raw_data_import[,3:length(colnames(raw_data_import))]

#Fix column names
colnames(raw_data)[which(colnames(raw_data)=="X184A1")]="184A1"
colnames(raw_data)[which(colnames(raw_data)=="X184B5")]="184B5"
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"
colnames(raw_data)[which(colnames(raw_data)=="X21MT1")]="21MT1"
colnames(raw_data)[which(colnames(raw_data)=="X21MT2")]="21MT2"
colnames(raw_data)[which(colnames(raw_data)=="X21NT")]="21NT"
colnames(raw_data)[which(colnames(raw_data)=="X21PT")]="21PT"

#Filter down to core set of cell lines - must have drug data and at least one other molecular profiling data type
core_cell_lines=c("600MPE","AU565","BT20","BT474","BT483","BT549","CAMA1","HCC1143","HCC1187","HCC1395","HCC1419","HCC1428","HCC1569","HCC1806","HCC1937","HCC1954","HCC202","HCC2185","HCC3153","HCC38","HCC70","HS578T","LY2","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","SKBR3","SUM1315MO2","SUM149PT","SUM159PT","SUM185PE","SUM225CWN","SUM44PE","SUM52PE","T47D","UACC812","UACC893","ZR751","ZR7530","ZR75B")
data=raw_data[,core_cell_lines]
cell_line_data=raw_cell_line_import[core_cell_lines,]

#Set row names
rownames(data)=paste(raw_feat_data[,"DataType"],raw_feat_data[,"ID"], sep="__")

#Categorical variables
cat_vars=which(raw_feat_data[,"DataType"]=='ES' | raw_feat_data[,"DataType"]=='MM' | raw_feat_data[,"DataType"]=='CNV')

#Use RandomForests to Impute missing values. Use subtype as the target so that cell lines of the same type will be assigned more similar values
target=as.factor(as.vector(cell_line_data[,"BCCLclassification"]))

#Transpose data so that predictor variables are columns and then set categorical variables to factors
tdata=as.data.frame(t(data))
tdata[,cat_vars]=sapply(tdata[,cat_vars],as.factor, simplify=FALSE)

data_imputed=rfImpute(x=tdata, y=target, ntree=300, iter=5)#Increase numbers of ntree and iter for final result?
data_imputed_clean=t(data_imputed[,2:length(colnames(data_imputed))])
data_imputed_clean=cbind(raw_feat_data,data_imputed_clean)

write.table(data_imputed_clean,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


