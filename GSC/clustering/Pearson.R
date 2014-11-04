data1 <- scan("/home/obig/tmp/value_tempfile",list(x=0,y=1))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x1 <- data1$x
x2 <- data1$y

r <- cor(x1, x2, method=c("pearson")) #calculates a Pearson correlation
p <- cor.test(x1, x2, method = "pearson", alternative = "greater")

sink("/home/obig/tmp/Pearson_tempfile") #sink directs output to specified file
print (r)
print (p)
sink() #returns output to R prompt