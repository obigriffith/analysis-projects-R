#datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/Affy"
datadir = "~/Projects/cisRED_coexpression/commonModuleAnalysis/cisred_2/TMM"
setwd(datadir)

#Postscript/pdf files
#ps_file="LogFreq_dist_rand_vs_norm.ps"
#pdf(file="CommonModulesAffy_vs_random.pdf", width=7, height=7)
#pdf(file="CommonModules_AB2_vs_AB3plus_vs_DN_vs_random_AFFY.pdf", width=7, height=7)
#pdf(file="CommonModules_AB2_vs_AB3plus_vs_DN_vs_random_TMM.pdf", width=7, height=7)

#Get Coexpression data for genes defined by modules or random
#Affy
##data_rand = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_rand = read.table("CommonModuleGeneClusters_allgenes.Affy_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB2 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_doubletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
##data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPlusPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
#data_norm_DN = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.Affy_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)

#TMM
data_rand = read.table("CommonModuleGeneClusters_allgenes.TMM_by_clustersize.random10000pairs_20.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_AB2 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_doubletPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_AB3 = read.table("AnnotationBased/CommonModuleGeneClusters_AnnotationBased_tripletPlusPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)
data_norm_DN = read.table("profileThresholded_250K/CommonModuleGeneClusters_profileThresholded_250K_tripletPatterns.TMM_by_clustersize.txt", header=F, quote="", sep="\t", comment.char="", as.is=1)


#Specify a maximum cluster size
maxclustersize=20
#maxclustersize=10
data_norm_AB3_undermax = data_norm_AB3[data_norm_AB3[,1]<=maxclustersize,]
data_norm_AB2_undermax = data_norm_AB2[data_norm_AB2[,1]<=maxclustersize,]
data_norm_DN_undermax = data_norm_DN[data_norm_DN[,1]<=maxclustersize,]

##To plot frequency of all Pearson for random vs normal (multiple methods can each be plotted)
#Get numbers of datapoints to calculate frequencies
num_coexp_norm_AB3 = length(data_norm_AB3_undermax[,2])
num_coexp_norm_AB2 = length(data_norm_AB2_undermax[,2])
num_coexp_norm_DN = length(data_norm_DN_undermax[,2])
num_coexp_rand = length(data_rand[,2])

#Sort coexpression scorres
coexp_norm_AB3_sorted=sort(data_norm_AB3_undermax[,2])
coexp_norm_AB2_sorted=sort(data_norm_AB2_undermax[,2])
coexp_norm_DN_sorted=sort(data_norm_DN_undermax[,2])
coexp_rand_sorted=sort(data_rand[,2])

#Calculate cumulative fraction of genepairs
fractions_norm_AB3=(1:num_coexp_norm_AB3)/num_coexp_norm_AB3
fractions_norm_AB2=(1:num_coexp_norm_AB2)/num_coexp_norm_AB2
fractions_norm_DN=(1:num_coexp_norm_DN)/num_coexp_norm_DN
fractions_rand=(1:num_coexp_rand)/num_coexp_rand


#Prepare to plot cum fractions for each on one plot
#plot(c(coexp_norm_AB3_sorted,coexp_norm_DN_sorted,coexp_rand_sorted), c(fractions_norm_AB3,fractions_norm_DN,fractions_rand), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Coexpression score", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
#plot(c(coexp_norm_AB2_sorted,coexp_norm_AB3_sorted,coexp_rand_sorted), c(fractions_norm_AB2,fractions_norm_AB3,fractions_rand), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Coexpression score", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
plot(c(coexp_norm_AB2_sorted,coexp_norm_AB3_sorted,coexp_norm_DN_sorted,coexp_rand_sorted), c(fractions_norm_AB2,fractions_norm_AB3,fractions_norm_DN,fractions_rand), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Coexpression score", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)
#plot(c(coexp_norm_AB2_sorted,coexp_norm_AB3_sorted,coexp_norm_DN_sorted,coexp_rand_sorted), c(fractions_norm_AB2,fractions_norm_AB3,fractions_norm_DN,fractions_rand), type="n", xlim=c(0,25), ylim=c(0,1), xlab="Coexpression score", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)

#Add points
points(coexp_norm_AB2_sorted, fractions_norm_AB2, pch=19, cex=0.25, col="orange")
points(coexp_norm_AB3_sorted, fractions_norm_AB3, pch=22, cex=0.25, col="darkgreen")
points(coexp_norm_DN_sorted, fractions_norm_DN, pch=23, cex=0.25, col="blue")
points(coexp_rand_sorted, fractions_rand, pch=25, cex=0.25, col="red")

#Plot a legend
#legend=c("Annotation Based (triplet)", "De Novo", "Random")
#legend=c("Annotation Based (doublet)", "Annotation Based (triplet)", "Random")
#legend=c("Annotation Based (doublet)", "Annotation Based (triplet)", "De Novo (triplet)", "Random")
legend=c("Annotation Based (doublet)", "Annotation Based (triplet+)", "De Novo (triplet)", "Random")

#legend(0.4,0.4, legend, pch=c(19,23,25), col=c("darkgreen","blue","red"), cex=1.5)
#legend(0.3,0.4, legend, pch=c(19,23,25), col=c("orange","darkgreen","red"), cex=1.5)
legend(0.3,0.4, legend, pch=c(20,15,18,17), col=c("orange","darkgreen","blue","red"), cex=1.25)

dev.off()

#Test for significance between distributions
#Problem is this tie warning an issue?
ks.test(x=coexp_norm_AB3_sorted, y=coexp_rand_sorted, alternative = c("two.sided"), exact = NULL)
ks.test(x=coexp_norm_DN_sorted, y=coexp_rand_sorted, alternative = c("two.sided"), exact = NULL)

