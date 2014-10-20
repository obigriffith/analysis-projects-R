#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/intellectual_property/RFRS/")

datafile="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/processing/processed_final2/test_survival/combined/ALL_gcrma.txt" #combined (standardCDF + customCDF)

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
header=colnames(data_import)

#Get predictor variables
top17opt_probes=c("204767_s_at","10682_at","201291_s_at","9133_at","1164_at","208079_s_at","23224_at","55435_at","23220_at","201461_s_at","202709_at","57122_at","23405_at","201483_s_at","29127_at","204416_x_at","10628_at")
top17opt_data=data_import[data_import[,1]%in%top17opt_probes,]
predictor_data=t(top17opt_data[,4:length(header)]) #Top20 optimized list
predictor_names=top17opt_data[,3] #gene symbol
colnames(predictor_data)=predictor_names

#Extract sample data for one patient to use as test case
#GSM36893 - predicted high risk (RFRS=0.810)

sample_predictor_data=as.data.frame(t(predictor_data["GSM36893.CEL",]))
rownames(sample_predictor_data)="GSM36893"

write.table(sample_predictor_data, file="patient_data.txt", col.names=NA, quote=FALSE, sep="\t")
