data1 <- scan("datafile.temp",list(x=0,y=1))

#Assign data from files to x0 and x1 from data1 file using the x and y labels specified in the list function above
x1 <- data1$x
x2 <- data1$y

r <- cor(x1, x2, method=c("pearson")) #calculates a Pearson correlation

sink("pearson_value.temp") #sink directs output to specified file
print (r)
sink() #returns output to R prompt

#use the following function if you also want a p-value associated with the pearson value.
#p <- cor.test(x1, x2, method = "pearson", alternative = "greater")
#print (p)

#options(echo = FALSE)