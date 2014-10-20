outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/filtered/"

setwd(outdir)

#stringent
pe_thresh1 = 0.2 #Minimum percent libraries "expressed"
cov_min1 = 2 #Minimum coefficient of variation (0.7 recommended?)
cov_max1 = 10 #Maximum cov

#weak
pe_thresh2 = 0.2 #Minimum percent libraries "expressed"
cov_min2 = 1 #Minimum coefficient of variation (0.7 recommended?)
cov_max2 = 10 #Maximum cov


#RNAseq data (Choose one feature type) 
#1. Gene
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_GeneExpression_v53.txt"
outfile1="breastRNAseq_genelevel_stringent.txt" #1A
outfile2="breastRNAseq_genelevel_weak.txt" #1B

#2. Transcript
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_TranscriptExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_TranscriptExpression_v53.txt"
outfile1="breastRNAseq_transcriptlevel_stringent.txt" #2A
outfile2="breastRNAseq_transcriptlevel_weak.txt" #2B

#3. Junction
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_JunctionExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_JunctionExpression_v53.txt"
outfile1="breastRNAseq_junctionlevel_stringent.txt" #3A
outfile2="breastRNAseq_junctionlevel_weak.txt" #3B

#4. Boundary
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_BoundaryExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_BoundaryExpression_v53.txt"
outfile1="breastRNAseq_boundarylevel_stringent.txt" #4A
outfile2="breastRNAseq_boundarylevel_weak.txt" #4B

#5. ExonRegion
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ExonRegionExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ExonRegionExpression_v53.txt"
outfile1="breastRNAseq_exonlevel_stringent.txt" #5A
outfile2="breastRNAseq_exonlevel_weak.txt" #5B

#6. ActiveIntron
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ActiveIntronRegionExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ActiveIntronRegionExpression_v53.txt"
outfile1="breastRNAseq_activeintronlevel_stringent.txt" #6A
outfile2="breastRNAseq_activeintronlevel_weak.txt" #6B

#7. ActiveIntergenic
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Matrix_ActiveIntergenicRegionExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/67libs/matrix/Expressed_ActiveIntergenicRegionExpression_v53.txt"
outfile1="breastRNAseq_activeintergeniclevel_stringent.txt" #6A
outfile2="breastRNAseq_activeintergeniclevel_weak.txt" #6B

#Import chosen data type
raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#Break data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:length(colnames(raw_data_import))]
raw_exp_status=raw_exp_status_import[,5:length(colnames(raw_data_import))]

#Cell lines to exclude from RNAseq data. Low qual: BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415; suspension line: DU4475; serum-free-media: SUM190PT; Normal; Other: M4A4, NM2C5; 
raw_data=subset(raw_data, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 
raw_exp_status=subset(raw_exp_status, select = -c(BT20, HCC1187, HCC1500, SUM159PT, SUM185PE, MDAMB157, MDAMB415, DU4475, SUM190PT, Normal, M4A4, NM2C5)) 

#Fix misnamed library
colnames(raw_data)[which(colnames(raw_data)=="MDAMB13v1")]="MDAMB134VI"
colnames(raw_data)[which(colnames(raw_data)=="SUM1315")]="SUM1315MO2"
colnames(raw_data)[which(colnames(raw_data)=="X184A1")]="184A1"
colnames(raw_data)[which(colnames(raw_data)=="X184B5")]="184B5"
colnames(raw_data)[which(colnames(raw_data)=="X600MPE")]="600MPE"
colnames(raw_data)[which(colnames(raw_data)=="X21MT1")]="21MT1"
colnames(raw_data)[which(colnames(raw_data)=="X21MT2")]="21MT2"
colnames(raw_data)[which(colnames(raw_data)=="X21NT")]="21NT"
colnames(raw_data)[which(colnames(raw_data)=="X21PT")]="21PT"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="MDAMB13v1")]="MDAMB134VI"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="SUM1315")]="SUM1315MO2"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X184A1")]="184A1"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X184B5")]="184B5"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X600MPE")]="600MPE"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X21MT1")]="21MT1"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X21MT2")]="21MT2"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X21NT")]="21NT"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="X21PT")]="21PT"

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

#apply stringent filters
pe_data=apply(raw_exp_status, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh1)
data_filt1=raw_data[passed_pe,]
feat_data_filt1=raw_feat_data[passed_pe,]

cov_data=apply(data_filt1, 1, cov_fun)
passed_cov = which(cov_data < cov_max1 & cov_data > cov_min1)
data_filt1=data_filt1[passed_cov,]
feat_data_filt1=feat_data_filt1[passed_cov,]

#apply weak filters
pe_data=apply(raw_exp_status, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh2)
data_filt2=raw_data[passed_pe,]
feat_data_filt2=raw_feat_data[passed_pe,]

cov_data=apply(data_filt2, 1, cov_fun)
passed_cov = which(cov_data < cov_max2 & cov_data > cov_min2)
data_filt2=data_filt2[passed_cov,]
feat_data_filt2=feat_data_filt2[passed_cov,]

#Write filtered results to file
write.table(cbind(feat_data_filt1,data_filt1), file=outfile1, sep="\t", row.names=FALSE)
write.table(cbind(feat_data_filt2,data_filt2), file=outfile2, sep="\t", row.names=FALSE)


