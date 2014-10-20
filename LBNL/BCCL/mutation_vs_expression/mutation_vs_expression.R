
outdir="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/mutation_vs_expression/"
setwd(outdir)

#Libraries that have both RNAseq and mutation data
RNAseq_mut=c("BT20","BT474","BT549","CAMA1","HCC1143","HCC1419","HCC1428","HCC1500","HCC1569","HCC1806","HCC202","HCC3153","HCC38","HCC70","HS578T","MCF10A","MCF10F","MCF12A","MCF7","MDAMB134VI","MDAMB157","MDAMB175VII","MDAMB231","MDAMB361","MDAMB453","SUM1315","SUM149PT","SUM159PT","T47D","UACC812","UACC893","ZR751","ZR7530","NM2C5")

#Load expression data, expression status, and mutation status
datafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Matrix_GeneExpression_v53.txt"
expressedfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/ALEXA/57libs/matrix/Expressed_GeneExpression_v53.txt"
mutdatafile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/mutations/ALL_filt_mutations.txt"
mutcountsfile="C:/Users/Obi/Documents/My Dropbox/Projects/BCCL/mutations/ALL_filt_mutations_counts.txt"

raw_data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_exp_status_import=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
raw_mut_import=read.table(mutdatafile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)
raw_mut_counts_import=read.table(mutcountsfile, header = TRUE, na.strings = "NA", sep="\t", row.names=1)

#Break RNAseq data into features info and expression values
raw_feat_data=raw_data_import[,1:4]
raw_data=raw_data_import[,5:62]
raw_exp_status=raw_exp_status_import[,5:62]

#Fix misnamed library
colnames(raw_data)[which(colnames(raw_data)=="MDAMB13v1")]="MDAMB134VI"
colnames(raw_exp_status)[which(colnames(raw_exp_status)=="MDAMB13v1")]="MDAMB134VI"

#Retrieve data for only libraries with both RNAseq and mutation data
data=raw_data[,RNAseq_mut]
feat_data=raw_feat_data
exp_status=raw_exp_status[,RNAseq_mut]
libs=colnames(data)
mut_data=raw_mut_import[,RNAseq_mut]
mut_counts=raw_mut_counts_import[,RNAseq_mut]


#Retrieve data for only genes of interest
#genes_interest=c("TP53","PIK3CA","ATM","PTEN","HRAS","CDH1","APC","BRAF","BRCA2","KRAS","AKT1","CDKN2A","RB1","ERBB2","EGFR","NF2","SMAD4","BRCA1","GATA3","NRAS","PIK3R1","JAK2")
#mut_data_interest=mut_data[genes_interest,]
#mut_counts_interest=mut_counts[genes_interest,]
#apply(mut_counts_interest,1,sum)

#retrieve data for only genes of interest that actually have at least one mutation in cell line that also has expression data
#genes_interest2=c("TP53","PIK3CA","ATM","HRAS","CDH1","PTEN","APC","BRAF","BRCA2","KRAS","CDKN2A","AKT1","RB1","ERBB2","EGFR","NF2")
#mut_data_interest2=mut_data[genes_interest2,]
#mut_counts_interest2=mut_counts[genes_interest2,]
#apply(mut_counts_interest2,1,sum)

#retrieve data for only genes of interest (Anneleen)
#genes_interest2=c("PTEN","BRCA1","BRCA2","PIK3CA")
#mut_data_interest2=mut_data[genes_interest2,]
#mut_counts_interest2=mut_counts[genes_interest2,]
#apply(mut_counts_interest2,1,sum)

#genes_interest2=c("ERBB2")
#mut_data_interest2=mut_data[genes_interest2,]
#mut_counts_interest2=mut_counts[genes_interest2,]
#apply(mut_counts_interest2,1,sum)

genes_interest2=c("NOTCH1")
mut_data_interest2=mut_data[genes_interest2,]
mut_counts_interest2=mut_counts[genes_interest2,]
apply(mut_counts_interest2,1,sum)


#Calculate statistics for mut
mut_totals=vector(mode="numeric", length=length(genes_interest2))
exp_means=vector(mode="numeric",length=length(genes_interest2))

for (i in 1:length(genes_interest2)){
 gene=genes_interest2[i]
 mut_data_gene=mut_data[gene,]
 exp_data_gene=log2(data[which(feat_data[,3]==gene),]+1)
 exp_status_gene=exp_status[which(feat_data[,3]==gene),]
 feat_data_gene=feat_data[which(feat_data[,3]==gene),]
 mut_totals[i]=sum(mut_data_gene)
 exp_means[i]=mean(as.numeric(exp_data_gene))
}

xlabels=paste(genes_interest2, " (", mut_totals, ")", sep="")



#Make summary figure showing for each gene, expression level, with points coded by mutation status
#Determine y max value
max_y=0
for (i in 1:length(genes_interest2)){
 gene=genes_interest2[i]
 exp_data_gene=log2(data[which(feat_data[,3]==gene),]+1)
 max=max(exp_data_gene)
 if(max>max_y){max_y=max}
}

#pdf(file="RNAseq_vs_mutation.pdf")
#pdf(file="RNAseq_vs_mutation_PTEN_BRCA_PI3K.pdf")
pdf(file="RNAseq_vs_mutation_NOTCH1.pdf")

plot(x=1, type="n", xlim=c(1,length(genes_interest2)), ylim=c(0,round(max_y+1.5)), xaxt="n", xlab="", ylab="Log2 Normalized RNAseq value", main="RNAseq expression versus mutation status (34 lines)")
axis(1, at=1:length(genes_interest2), labels=xlabels, las=2, cex.axis=0.8)
legend("topleft", legend=c("Exp, Mt", "Exp, Wt", "Not Exp, Mt","Not Exp, Wt" ), pch=c(19,4,19,4), col=c("red","red","blue","blue"))

for (i in 1:length(genes_interest2)){
 gene=genes_interest2[i]
 mut_data_gene=mut_data[gene,]
 exp_data_gene=log2(data[which(feat_data[,3]==gene),]+1)
 exp_status_gene=exp_status[which(feat_data[,3]==gene),]
 feat_data_gene=feat_data[which(feat_data[,3]==gene),]

 #Set colors and symbols to distinguish wt from mut
 mut_status_symbols=as.numeric(mut_data_gene)
 mut_status_symbols[mut_status_symbols==1]=19
 mut_status_symbols[mut_status_symbols==0]=4
 exp_status_colors=as.numeric(exp_status_gene)
 exp_status_colors[exp_status_colors==0]="blue"
 exp_status_colors[exp_status_colors==1]="red"

 points(x=rep(i,length(exp_data_gene)), y=exp_data_gene, pch=mut_status_symbols,col=exp_status_colors)
}
dev.off()

