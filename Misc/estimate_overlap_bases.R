
estimate_loss=function(n,mean,sd){
	data=rnorm(n=n,mean=mean,sd=sd)
	data_round=round(data)
	data_round_lt200=data_round[data_round<200]
	data_round_lt200_min50=data_round_lt200
	data_round_lt200_min50[data_round_lt200_min50<50]=50
	lost_bases=data_round_lt200_min50-200
	perc_lost=((sum(lost_bases)*-1)/(1000000*200))*100
	return(perc_lost)
}
n=1000000
mean=370
sd=76

estimate_loss(n,mean,sd)
