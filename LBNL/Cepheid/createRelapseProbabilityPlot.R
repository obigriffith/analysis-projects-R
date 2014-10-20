
#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/finaltop17/")

#datafile="Cepheid_CasePredictions_combined_downsamp.txt" #Single down-sampled iteration
datafile="Cepheid_CasePredictions_combined.txt" #All data, will have to down-sample

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")

#Extract data for plotting
case_predictions=data_import[,c("X10yr_relapse","Relapse")]

#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Keep all NoRelapses
NoRelapseCases=which(is.na(case_predictions[,"X10yr_relapse"]) | case_predictions[,"X10yr_relapse"]==0)
RelapseCases=which(case_predictions[,"X10yr_relapse"]==1)

#Downsample Relapse cases so that they represent only 15% of total cases:  x / (429+x) = 0.15 [solving for x, ~76]
#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
N=100

#Break 10yr relapse frequency into bins (by Relapse score) for plotting
max_score=max(case_predictions[,"Relapse"])
min_score=min(case_predictions[,"Relapse"])

n_bins=50
#bin_size=1/n_bins #relapse score is value from 0 to 1
bin_size=(max_score-min_score)/n_bins #relapse score is value from min to max

#breaks=seq(0,1,by=bin_size)
breaks=seq(min_score,max_score,by=bin_size)

#Create matrix to store results in
downsampledata=matrix(data=NA, nrow=N, ncol=length(breaks)-1)
colnames(downsampledata)=breaks[2:length(breaks)]

for (i in 1:N){
 random_RelapseCases=sample(x=RelapseCases, size=76, replace = FALSE, prob = NULL)
 case_predictions_down=case_predictions[c(NoRelapseCases,random_RelapseCases),]
 for (j in 1:length(breaks)){
  if(j<length(breaks)){
   #Get data for bin
   bin_relapses=case_predictions_down[which(case_predictions_down[,"Relapse"]>breaks[j] & case_predictions_down[,"Relapse"]<=breaks[j+1]),"X10yr_relapse"]
   bin_relapse_fraction=sum(bin_relapses,na.rm=TRUE)/length(bin_relapses)
   downsampledata[i,j]=bin_relapse_fraction
  }
 }
}

#Calculate means for N iterations
bin_relapse_fraction_down_means=apply(downsampledata, 2, mean)

#Convert relapse fractions to percentages
bin_relapse_fraction_down_means_perc=bin_relapse_fraction_down_means*100

#plot mean relapse fractions by bin, fit a line, and calculate confidence intervals
#First calculate fit lien
df=data.frame(mean_relapse_fraction=bin_relapse_fraction_down_means_perc, breaks=as.numeric(names(bin_relapse_fraction_down_means_perc)))
fit=loess(df$mean_relapse_fraction ~ df$breaks, span=0.9, degree=1)
pred=predict(fit, df$breaks, se=TRUE)

#Create plot
plot(df$mean_relapse_fraction ~ df$breaks, ylab="Likelihood of Relapse at 10 years (%)", xlab="Random Forests Relapse Score (RFRS)", type="n", ylim=c(0,50), xlim=c(0,1), xaxs = "i", las=1, cex.lab=1.2, cex.axis=1.2)

#Extrapolate more points to make line smoother
xl=seq(min(df$breaks),max(df$breaks), (max(df$breaks) - min(df$breaks))/1000)
xl.rev=sort(xl, decreasing=TRUE)
pred=predict(fit,xl,se=TRUE)
lines(x=xl, y=pred$fit, lty="solid", col="darkred", lwd=3)
lines(x=xl, y=pred$fit-1.96*pred$se.fit, lty="dashed", col="blue", lwd=1)
lines(x=xl, y=pred$fit+1.96*pred$se.fit, lty="dashed", col="blue", lwd=1)

###create polygon bounds
y.polygon=c((pred$fit+1.96*pred$se.fit)[1:length(xl)], (pred$fit-1.96*pred$se.fit)[length(xl):1])
x.polygon=c(xl, xl.rev)
polygon(x.polygon, y.polygon, col="#00009933", border=NA)

#Add vertical lines breaking plot into low, intermediate, and high risk group
abline(v=0.333)
abline(v=0.606)

#Add text indicating risk groups
text(x=0.17, y=45, labels="Low", pos=3, cex = 1.2)
text(x=0.47, y=45, labels="Intermediate", pos=3, cex = 1.2)
text(x=0.8, y=45, labels="High", pos=3, cex = 1.2)

#Add grid lines
grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

#Add tick marks to show distribution of Relapse scores for entire patient population
xsegs=case_predictions[,"Relapse"]
y0segs=rep(-2.5,length(xsegs))
y1segs=rep(0,length(xsegs))
segments(x0=xsegs, y0=y0segs, x1=xsegs, y1=y1segs)
RelapseProbabilityPlot=recordPlot() 
dev.off()

pdf(file="RFRS_vs_EstRelapseProbability.pdf")
replayPlot(RelapseProbabilityPlot)
dev.off()

save(RelapseProbabilityPlot,file="RelapseProbabilityPlot.Rdata")
save(fit,file="RelapseProbabilityFit.Rdata")

