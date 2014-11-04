#postscript("/home/obig/clustering/R_analysis/freq_dists/Affy_cum_freq_dist.ps", pointsize=1, horizontal=TRUE, width=5, height=4)
x <- scan("/home/obig/clustering/R_analysis/temp/affy_only_gt100_subset1_values.txt")
#data1 <- scan("/home/obig/clustering/R_analysis/temp/affy_only_gt100_subset1.txt",list(x=2))
#data1 <- scan("/home/obig/clustering/R_analysis/temp/affy_only_gt100_subset2.txt",list(x=2))
#data1 <- scan("/home/sage/clustering/matrix/affy_micro_sage_comp/040407/results/for_coexpression_db/affy_only_gt100.txt",list(x=2))
#x <- data1$x

#y <- scan("/home/grobertson/R/KolmogorovSmirnov/20040825/ENSG00000107796_FAKE.scores")
#test <- ks.test(x, y)
#pval <- test$p.value
# distribution function plot

edfx <- ( 1:length(x) )/length(x)

#edfy <- ( 1:length(y) )/length(y)
#par(cex=1.5)
#maxx <- max( max(x),max(y) )
#minx <- min( min(x),min(y) )
#xrange <- c(minx,maxx)
xrange <- c(-1,1)
#yrange <- c(0,0.4)
par (cex=2)
px <- plot(sort(x),edfx, type="s",
  xlab="Pearson Correlation",
  ylab="Distribution function", lwd=1.5, 
  xlim=xrange
)

#  ylim=yrange
#)


#text( min(sx)+0.05*range(sx), max(edfx)-0.05*range(edfx),
#  "Kolmogorov-Smirnov p-value =", pos=4)
#text( min(x)+0.05*range(x), min(edfx)+0.05*range(edfx), "ensg", pos=4)
#text( min(y)+0.05*range(y), min(edfy)+0.1*range(edfy), "fake", pos=4 )
#par(new=T)
#lines(sort(y),edfy,type="s", lty=8, lwd=1.5, col="gray")
#grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted", lwd = NULL, equilogs = TRUE)

#dev.print(device=pdf,file="/home/obig/clustering/R_analysis/freq_dists/Affy_cum_freq_dist.pdf")
dev.copy(device=png, file="/home/obig/clustering/R_analysis/freq_dists/Affy_cum_freq_dist.png", width=600, height=600)
dev.off()


# Does x come from a shifted gamma distribution with shape 3 and scale 2?
#ks.test(rx+2, "pgamma", 3, 2) # two-sided
#ks.test(rx+2, "pgamma", 3, 2, alternative = "gr")
# shapiro test for normality
#shapiro.test(rx)
#shapiro.test(ry)
