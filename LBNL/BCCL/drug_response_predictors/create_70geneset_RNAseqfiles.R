outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/70geneset/"
setwd(outdir)
outfile="breastRNAseq_all_70geneset.txt"

#Parameters
#stringent
pe_thresh1 = 0.2 #Minimum percent libraries "expressed"
cov_min1 = 2 #Minimum coefficient of variation (0.7 recommended?)
cov_max1 = 10 #Maximum cov

#Define a percent expressed function
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
  return(cov)
}

#Get list of ~70 genes from RPPA data
genelistfile="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RPPA/filtered/RPPA_GeneSymbol_mapping.txt"
genelist_import=read.table(genelistfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2))
gene_list=unique(genelist_import[,"Gene_Symbols"])

#Note, MTOR is called FRAP1 in the version of Ensembl used in RNAseq Alexa-seq analysis
gene_list[gene_list=="MTOR"]="FRAP1"

#Get corresponding ENSG ids for gene symbols
RNA_annotation_file="C:/Users/Obi/Documents/My Dropbox/drug_predictors/RNAseq/hs_53_36o_ENSG2Symbol.txt"
RNA_anno_import=read.table(RNA_annotation_file, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:2))

#Look up Symbol for ENSG ID and add to feature data
Symbol2ENSG=matrix(data=NA, nrow=length(gene_list), ncol=2)
for (i in 1:length(gene_list)){
  mapping=RNA_anno_import[RNA_anno_import[,"Gene_Name"]==gene_list[i],]
  Symbol2ENSG[i,1]=as.character(mapping[1])
  Symbol2ENSG[i,2]=as.character(mapping[2])
}

ENSG_list=Symbol2ENSG[,1]

#RNAseq data
#1. Gene
genedatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
#geneexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"
genedata_import=read.table(genedatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
genedata_set=genedata_import[genedata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]

#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
genedata_set=subset(genedata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 

#Create a common header for all datatypes to allow final rbind
header=colnames(genedata_set)
header[1]="ID"
colnames(genedata_set)=header

#Change FRAP1 back to MTOR for consistency with other datasets
genedata_set[genedata_set[,"Seq_Name"]=="FRAP1","Seq_Name"]="MTOR"


#Filter non-gene features to reduce numbers of pointless features (keep all genes)
#2. Transcript
transcriptdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_TranscriptExpression_v53.txt"
transcriptexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_TranscriptExpression_v53.txt"
transcriptdata_import=read.table(transcriptdatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
transcriptdata_set=transcriptdata_import[transcriptdata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
transcriptexpressed_import=read.table(transcriptexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
transcriptexpressed_set=transcriptexpressed_import[transcriptexpressed_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
transcriptdata_set=subset(transcriptdata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
transcriptexpressed_set=subset(transcriptexpressed_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
colnames(transcriptdata_set)=header

transcript_feat_data=transcriptdata_set[,1:4]
transcriptdata=transcriptdata_set[,5:length(header)]
transcriptexpressed=transcriptexpressed_set[,5:length(header)]

#apply stringent filters
pe_data=apply(transcriptexpressed, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=transcriptdata[passed_pe,]
feat_data_filt1=transcript_feat_data[passed_pe,]
cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]
transcriptdata_filt=cbind(feat_data_filt1,data_filt1)

#3. Junction
junctiondatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_JunctionExpression_v53.txt"
junctionexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_JunctionExpression_v53.txt"
junctiondata_import=read.table(junctiondatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
junctiondata_set=junctiondata_import[junctiondata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
junctionexpressed_import=read.table(junctionexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
junctionexpressed_set=junctionexpressed_import[junctionexpressed_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
junctiondata_set=subset(junctiondata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
junctionexpressed_set=subset(junctionexpressed_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
colnames(junctiondata_set)=header

junction_feat_data=junctiondata_set[,1:4]
junctiondata=junctiondata_set[,5:length(header)]
junctionexpressed=junctionexpressed_set[,5:length(header)]

#apply stringent filters
pe_data=apply(junctionexpressed, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=junctiondata[passed_pe,]
feat_data_filt1=junction_feat_data[passed_pe,]
cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]
junctiondata_filt=cbind(feat_data_filt1,data_filt1)

#4. Boundary
boundarydatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_BoundaryExpression_v53.txt"
boundaryexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_BoundaryExpression_v53.txt"
boundarydata_import=read.table(boundarydatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
boundarydata_set=boundarydata_import[boundarydata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
boundaryexpressed_import=read.table(boundaryexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
boundaryexpressed_set=boundaryexpressed_import[boundaryexpressed_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
boundarydata_set=subset(boundarydata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
boundaryexpressed_set=subset(boundaryexpressed_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
colnames(boundarydata_set)=header

boundary_feat_data=boundarydata_set[,1:4]
boundarydata=boundarydata_set[,5:length(header)]
boundaryexpressed=boundaryexpressed_set[,5:length(header)]

#apply stringent filters
pe_data=apply(boundaryexpressed, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=boundarydata[passed_pe,]
feat_data_filt1=boundary_feat_data[passed_pe,]
cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]
boundarydata_filt=cbind(feat_data_filt1,data_filt1)

#5. ExonRegion
exondatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ExonRegionExpression_v53.txt"
exonexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ExonRegionExpression_v53.txt"
exondata_import=read.table(exondatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
exondata_set=exondata_import[exondata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
exonexpressed_import=read.table(exonexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
exonexpressed_set=exonexpressed_import[exonexpressed_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
exondata_set=subset(exondata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
exonexpressed_set=subset(exonexpressed_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
colnames(exondata_set)=header

exon_feat_data=exondata_set[,1:4]
exondata=exondata_set[,5:length(header)]
exonexpressed=exonexpressed_set[,5:length(header)]

#apply stringent filters
pe_data=apply(exonexpressed, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=exondata[passed_pe,]
feat_data_filt1=exon_feat_data[passed_pe,]
cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]
exondata_filt=cbind(feat_data_filt1,data_filt1)

#6. ActiveIntron
introndatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ActiveIntronRegionExpression_v53.txt"
intronexpressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ActiveIntronRegionExpression_v53.txt"
introndata_import=read.table(introndatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
introndata_set=introndata_import[introndata_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
intronexpressed_import=read.table(intronexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
intronexpressed_set=intronexpressed_import[intronexpressed_import[,"EnsEMBL_Gene_ID"]%in%ENSG_list,]
#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
introndata_set=subset(introndata_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
intronexpressed_set=subset(intronexpressed_set, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
colnames(introndata_set)=header

intron_feat_data=introndata_set[,1:4]
introndata=introndata_set[,5:length(header)]
intronexpressed=intronexpressed_set[,5:length(header)]

#apply stringent filters
pe_data=apply(intronexpressed, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=introndata[passed_pe,]
feat_data_filt1=intron_feat_data[passed_pe,]
cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]
introndata_filt=cbind(feat_data_filt1,data_filt1)

#7. ActiveIntergenic  #Skip Intergenics because can't be simply mapped to gene symbols in the set
#datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ActiveIntergenicRegionExpression_v53.txt"
#expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ActiveIntergenicRegionExpression_v53.txt"
#outfile1="breastRNAseq_activeintergeniclevel_stringent.txt" #6A
#outfile2="breastRNAseq_activeintergeniclevel_weak.txt" #6B

#Combine all features into a single dataset
alldata_set=rbind(genedata_set, transcriptdata_filt, exondata_filt, junctiondata_filt, boundarydata_filt, introndata_filt)

#Fix misnamed library
colnames(alldata_set)[which(colnames(alldata_set)=="MDAMB13v1")]="MDAMB134VI"
colnames(alldata_set)[which(colnames(alldata_set)=="X184A1")]="184A1"
colnames(alldata_set)[which(colnames(alldata_set)=="X184B5")]="184B5"
colnames(alldata_set)[which(colnames(alldata_set)=="X600MPE")]="600MPE"
colnames(alldata_set)[which(colnames(alldata_set)=="X21MT1")]="21MT1"
colnames(alldata_set)[which(colnames(alldata_set)=="X21MT2")]="21MT2"
colnames(alldata_set)[which(colnames(alldata_set)=="X21NT")]="21NT"
colnames(alldata_set)[which(colnames(alldata_set)=="X21PT")]="21PT"

write.table(alldata_set, file=outfile, sep="\t", row.names=FALSE)


