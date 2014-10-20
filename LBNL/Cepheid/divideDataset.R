clindatafile="C:/Users/Obi/Documents/Projects/Cepheid/clinical_data/filtered_combined_data_anno.final.txt"
outdir="C:/Users/Obi/Documents/Projects/Cepheid/clinical_data"
clindatafile_random_train="filtered_combined_data_anno.final.train.txt"
clindatafile_random_test="filtered_combined_data_anno.final.test.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#Change to output dir
setwd(outdir)
header=colnames(clin_data_import)

studies=levels(clin_data_import[,"Study"])

#Initialize empty vectors to store index values for training/test sets
training=c()
test=c()

#Go through each study and divide into training and test datasets
#Balance for X10yr_relapse
#Studies with NA for X10yr_relapse should be included in training data (can be used in Survival analyis)
for (i in 1:length(studies)){
  data_sub=clin_data_import[clin_data_import[,"Study"]==studies[i],]

  #Split data further into no10yrFU, 10yrFUrelapse, 10yrFUnorelapse
  data_sub_no10yrFU=data_sub[which(is.na(data_sub[,"X10yr_relapse"])),]
  data_sub_10yrFUrelapse=data_sub[which(data_sub[,"X10yr_relapse"]==1),]
  data_sub_10yrFUnorelapse=data_sub[which(data_sub[,"X10yr_relapse"]==0),]

  #Add no10yrFU to training data
  training=c(training,as.numeric(rownames(data_sub_no10yrFU)))

  #break 10yrFUrelapse/10yrFUnorelapse into 2/3 training 1/3 test
  #First, count number of relapses and no relapses for study
  num_10yrFUrelapse=length(rownames(data_sub_10yrFUrelapse))
  num_10yrFUnorelapse=length(rownames(data_sub_10yrFUnorelapse))

  #Determine what 1/3 of this number would be to put aside for testing
  num_test_10yrFUrelapse=round(num_10yrFUrelapse*1/3)
  num_test_10yrFUnorelapse=round(num_10yrFUnorelapse*1/3)

  #Randomly sample from 1 to total number of samples. Then take first 1/3 as test, remaining 2/3 as training
  random_10yrFUrelapse=sample(x=num_10yrFUrelapse, size=num_10yrFUrelapse, replace = FALSE, prob = NULL)
  random_10yrFUnorelapse=sample(x=num_10yrFUnorelapse, size=num_10yrFUnorelapse, replace = FALSE, prob = NULL)

  #Reorder data subsets according to randomization
  test=c(test, as.numeric(rownames(data_sub_10yrFUrelapse[random_10yrFUrelapse[1:num_test_10yrFUrelapse],])))
  training=c(training, as.numeric(rownames(data_sub_10yrFUrelapse[random_10yrFUrelapse[(num_test_10yrFUrelapse+1):num_10yrFUrelapse],])))
  test=c(test, as.numeric(rownames(data_sub_10yrFUnorelapse[random_10yrFUnorelapse[1:num_test_10yrFUnorelapse],])))
  training=c(training, as.numeric(rownames(data_sub_10yrFUnorelapse[random_10yrFUnorelapse[(num_test_10yrFUnorelapse+1):num_10yrFUnorelapse],])))
}

write.table(clin_data_import[test,], file=clindatafile_random_test, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(clin_data_import[training,], file=clindatafile_random_train, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



