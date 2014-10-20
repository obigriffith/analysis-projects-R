#Transcript-level fpkm data from cufflinks
cuffdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/formattedFiles/Her2_fpkm_v2.txt"

outdir="C:/Users/Obi/Documents/My Dropbox/Projects/SU2C/HER2/TophatCufflinks/formattedFiles/"
setwd(outdir)

#Import data
cuffdata_import=read.table(cuffdatafile, header = TRUE, na.strings = "NA", as.is=c(1:2))
cuffexpdata=cuffdata_import[,3:length(colnames(cuffdata_import))]

#Merge/average/add at gene symbol level
cuffexpdata_genesum=aggregate.data.frame(cuffexpdata, by=list(cuffdata_import[,"GeneSymbol"]), sum)
colnames(cuffexpdata_genesum)[1]="GeneSymbol"

write.table(cuffexpdata_genesum, file="Her2_fpkm_v2_genelevel.txt", quote=FALSE, sep="\t", row.names=FALSE)