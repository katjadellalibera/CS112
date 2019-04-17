#######################
#   John Henderson    #
#    UC-Berkeley      #            
#       &             #
#   Sara Chatfield    #
#    UC-Berkeley      #       
#                     # 
#    Who Matches:     #
#   Replication Code  #
#                     #
#   Released          #
#     Version 1.0     #
#                     #
#   Nov 17, 2010      #
#######################


# USAGE
#########################################  
# This is replication code for Henderson, John, and Sara Chatfield. 2011. "Who Matches? Propensity Scores and Bias in the Causal Effects
# of Education on Participation," Journal of Politics.  The software, code, and data may be used, distributed, and modified freely 
# with proper attribution to both Henderson and Chatfield (2011) and Kam and Palmer (2008). 


# DESCRPITION
#########################################  
# This code produces the 766642 x 81 indicator matrix that is the engine behind the pscore iterations analysis in the paper.  
# The sampling method samples uniformly within strata defined at each n discrete factorial level of 81 choose n covariates. In 
# other words, for each n = 1, 2, ... 81, the sasmpling method samples up to 10,000 permutations of fulle set of 81 choose n 
# permutations.  This file is included mainly to clarify the sampling procedure described in prose detail in "Who Matches?".

# 1. Indicator

#81 choose N Sampling, for N = 3, 4, ..., 78.  

set.seed(1005)

items=list()
indicators=matrix()
for(k in 1:76)
{
	items[[k]]=matrix(0,nrow=10000,ncol=81)
	
	
	for(i in 1:10000)
	{
		items[[k]][i,]=sample(c(rbind(matrix(1,k+2),matrix(0,81-k-2))),81,replace=F)
	}
}

#81 choose 1

ind=matrix(0,nrow=81,ncol=81)
for(i in 1:81){
	ind[i,i]=1
}

indicator.1=ind

#81 choose 2

ind.list=list()
j=c(80:1)
vec.store=list()
vec=matrix()
out=list()
outs=matrix()

for(i in 1:80){
	
	temp=matrix(0,nrow=j[i],ncol=j[i])
	
	for(k in 1:j[i]){
		temp[k,k]=1
	}
	ind.list[[i]]=temp		
	
	for(l in 1:80){
		vec=matrix(0,nrow=j[l],ncol=80-j[l])
		vec=cbind(vec,matrix(1,nrow=j[l],ncol=1))
		vec.store[[l]]=vec
	}
	out[[i]]=cbind(vec.store[[i]],ind.list[[i]])
}	

outs=out[[1]]
for(g in 2:length(out)){
	outs=rbind(outs,out[[g]])
}

indicator.2=outs

#81 choose 79

ind.list=list()
j=c(80:1)
vec.store=list()
vec=matrix()
out=list()
outs=matrix()

for(i in 1:80){
	
	temp=matrix(1,nrow=j[i],ncol=j[i])
	
	for(k in 1:j[i]){
		temp[k,k]=0
	}
	ind.list[[i]]=temp		
	
	for(l in 1:80){
		vec=matrix(1,nrow=j[l],ncol=80-j[l])
		vec=cbind(vec,matrix(0,nrow=j[l],ncol=1))
		vec.store[[l]]=vec
	}
	out[[i]]=cbind(vec.store[[i]],ind.list[[i]])
}	

outs=out[[1]]
for(g in 2:length(out)){
	outs=rbind(outs,out[[g]])
}

indicator.3=outs


#81 choose 80

ind=matrix(1,nrow=81,ncol=81)
for(i in 1:81){
	ind[i,i]=0
}

indicator.4=ind

# Cant malloc in this loop!
#indicator=rbind(indicator.1,indicator.2)
#for(i in 1:76){
#	indicator=rbind(indicator,items[[i]])
#}
#indicator=rbind(indicator,indicator.3,indicator.4)

indicator=rbind(indicator.1,indicator.2,items[[1]],items[[2]],items[[3]],items[[4]],items[[5]],items[[6]],items[[7]],items[[8]],items[[9]],items[[10]],items[[11]],items[[12]],items[[13]],items[[14]],items[[15]],items[[16]],items[[17]],items[[18]],items[[19]],items[[20]],items[[21]],
	items[[22]],items[[23]],items[[24]],items[[25]],items[[26]],items[[27]],items[[28]],items[[29]],items[[30]],items[[31]],items[[32]],items[[33]],items[[34]],items[[35]],items[[36]],items[[37]],items[[38]],items[[39]],items[[40]],items[[41]],items[[42]],items[[43]],items[[44]],
	items[[45]],items[[46]],items[[47]],items[[48]],items[[49]],items[[50]],items[[51]],items[[52]],items[[53]],items[[54]],items[[55]],items[[56]],items[[57]],items[[58]],items[[59]],items[[60]],items[[61]],items[[62]],items[[63]],items[[64]],items[[65]],items[[66]],items[[67]],
	items[[68]],items[[69]],items[[70]],items[[71]],items[[72]],items[[73]],items[[74]],items[[75]],items[[76]],indicator.3,indicator.4)

# END