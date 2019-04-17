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
# The following code replicates the genetic matching analysis conducted in 'Who Matches'.  In training our genetic matching
# runs, we made a series of customizations aimed at optimizing balance in covariates. The first of these is the use of a different
# loss function than that which is used in the standard GenMatch package.  We optimize the proportion of all covariates that are 
# balanced across treatment and control groups at the p > .1 level. Note that these are all t-tests since all covariates are 
# dichotomized following Kam and Palmer (2008).  This loss function is 'lower.bound' and is described in greater detail in the 
# funcs.R file.      
#
# The second unique specification is that we use two different propensity scores as starting points in our genetic matching 
# runs.  A typical way to use propensity scores in a genetic matching run is to create a matrix X that includes the covariates 
# that GenMatch will match on, and attach the p-score to this matrix so that GenMatch will search for matching weights for the 
# X covariates plus a p-score P.  We simply attach P1 and P2 to X and match on this matrix.  To speed up the GenMatch runs, we 
# also orthogonalize X by regressing each covariate vector on P1 and then on P2, extracting the residuals to form the final matching 
# matrix.  Though this isnt necessary, it greatly speeds up the optimization procedure since it ensures that GenMatch only matches 
# on information in X that is orthogonal to P1 and P2. 
#
# Three final things: First, we check balance and match on a matrix object 'out', which is produced by dichotomizing an expanded set of covariates 
# that include the original 81 used in Kam and Palmer (2008) plus an additional 28 covariates.  Dichotomizing these covariates 
# produces a matrix of 1254 x 355.  Second, the results below represent the final product from a very long series of GenMatch 
# trainings.  The GenMatch weights objects 'starts' are the weights supplied by GenMatch from our final run, and thus are used here 
# to recover the final ATT and ATC estimates.  To replicate the actual GenMatch runs we conducted, you will need to create a new 
# 'starts' vector that has 100 in the first two elements and 0s in the remaining 355, and increase the 'pop.size' and 'max.generations' 
# to sufficiently large values to get improvements.  
# 
# Finally, this code sources both the funcs.R file and the objects.R file.  The former contains the relevant functions used here to 
# conduct our analysis, and the latter contains the code that transforms the original dataset into the various matrix objects used
# in GenMatch.  The indicator.Rdata contains the 'indicator' object that is the matrix used in the 766642 permutations of the 81 
# covariates.  Here it is modified to just include the best two ATT and two ATC propensity scores.
            

library(MASS)
library(Matching)
#load("WhoMatches.Rdata")
source('funcs.R')
source('objects.R')   
load('indicator.Rdata')     

floors=NULL
floor.values=NULL
reweighted = -1
attach(data_rep)
attach(factors)


# BEST GENMATCH

#ATT: prob 1 to 3

matTemp=as.matrix(original)%*%diag(indicator[1,])   #363708
temp=list()
for(i in 1:81){
	if(indicator[1,i]==1){
		temp[[i]]=factor(matTemp[,i])
	}
	else{
		temp[[i]]=matTemp[,i]
	}
}

pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+
	temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+
	temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+
	temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
etahat1=pscore$fitted.values


matTemp=as.matrix(original)%*%diag(indicator[2,])    #388106
temp=list()
for(i in 1:81){
	if(indicator[2,i]==1){
		temp[[i]]=factor(matTemp[,i])
	}
	else{
		temp[[i]]=matTemp[,i]
	}
}

pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+
	temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+
	temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+
	temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
etahat2=pscore$fitted.values

Xmat=cbind(etahat1,out)
for(i in 2:ncol(Xmat)){
	Xmat[,i]=lm(Xmat[,i]~etahat1)$residuals
}

Xmat=cbind(etahat2,Xmat)
for(i in 3:ncol(Xmat)){
	Xmat[,i]=lm(Xmat[,i]~etahat2)$residuals
}


starts=c(99.99925,100.0053,0.005715251,0.006408247,2.980232e-08,2.980232e-08,0.006262702,0.008668285,2.980232e-08,0.002887763,2.980232e-08,0.0009998784,0.002455857,0.003635949,2.980232e-08,2.980232e-08,2.980232e-08,0.00787601,0.006441565,
	0.006882892,0.004390646,0.001944808,0.005531695,2.980232e-08,0.008536832,0.006503187,2.980232e-08,0.009334766,0.00209233,2.980232e-08,2.980232e-08,2.980232e-08,0.009141416,0.006710892,2.980232e-08,0.00375902,2.980232e-08,2.980232e-08,
	0.002394659,0.002549385,2.980232e-08,0.001481797,0.008610365,2.980232e-08,2.980232e-08,2.980232e-08,0.009984436,0.006574805,0.002360336,0.00822486,2.980232e-08,2.980232e-08,0.004001202,0.007789226,2.980232e-08,0.008463702,0.006892591,
	2.980232e-08,0.003049413,2.980232e-08,0.001608269,0.004021724,0.009204788,0.001364288,0.001575996,2.980232e-08,0.008038192,2.980232e-08,0.007762872,2.980232e-08,0.005283075,0.003295967,2.980232e-08,2.980232e-08,0.004179459,173.6054,
	0.002824357,0.001318068,0.008916264,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.007247979,2.980232e-08,0.006996465,2.980232e-08,0.009502948,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.007643543,2.980232e-08,
	2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.004981017,0.000132719,0.00259975,0.007494614,855.7601,0.003140685,2.980232e-08,2.980232e-08,2.980232e-08,0.005202841,2.980232e-08,0.009273753,2.980232e-08,0.00141859,2.980232e-08,
	0.002784641,0.008928815,2.980232e-08,2.980232e-08,0.005609324,2.980232e-08,0.007712388,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.001323451,0.004780747,0.001246533,2.980232e-08,0.008068434,2.980232e-08,2.980232e-08,0.001275443,
	0.003073698,0.003277063,2.980232e-08,2.980232e-08,2.980232e-08,0.004064664,0.008520083,2.980232e-08,0.006505253,0.005767007,0.00651125,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.004590578,0.005835394,2.980232e-08,2.980232e-08,
	2.980232e-08,0.003976639,2.980232e-08,2.980232e-08,2.980232e-08,0.00442357,2.980232e-08,0.0001809476,2.980232e-08,2.980232e-08,0.00434598,0.004373333,2.980232e-08,0.001005735,0.009737443,2.980232e-08,2.980232e-08,0.002011738,2.980232e-08,
	0.006726748,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.0055867,42.69455,2.980232e-08,0.000141961,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.001291604,0.007333195,2.980232e-08,
	0.008610143,0.0007476773,2.980232e-08,2.980232e-08,0.004024254,0.00452202,0.001858119,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.009335303,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,
	0.007920756,2.980232e-08,0.00153351,0.006578747,2.980232e-08,2.980232e-08,2.980232e-08,0.003463541,2.980232e-08,2.980232e-08,2.980232e-08,603.442,0.001144914,0.007914557,2.980232e-08,2.980232e-08,0.003372544,0.007569384,0.004901839,
	2.980232e-08,2.980232e-08,0.0008925979,2.980232e-08,0.005461286,0.008944928,0.0009883558,0.003052224,0.003033079,2.980232e-08,0.008314635,0.0003773529,2.980232e-08,2.980232e-08,0.009716536,288.4406,2.980232e-08,2.980232e-08,2.980232e-08,
	2.980232e-08,0.008213271,69.9943,2.980232e-08,0.002875351,0.00960413,0.007875248,0.005444037,2.980232e-08,0.003010958,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.006345055,2.980232e-08,0.005951712,0.0086756,2.980232e-08,2.980232e-08,
	0.0003141093,2.980232e-08,2.980232e-08,0.009148151,0.004543724,2.980232e-08,2.980232e-08,0.003030723,0.006520301,0.008936021,0.002618584,0.000648524,2.980232e-08,0.009145669,0.002337948,2.980232e-08,0.008956889,0.009425506,0.004757528,
	0.006359484,2.980232e-08,0.0002088765,0.007746522,2.980232e-08,2.980232e-08,2.980232e-08,392.4657,367.7674,2.980232e-08,0.001256075,2.980232e-08,2.980232e-08,0.009527099,0.009801216,0.002811012,0.006189428,0.00491197,2.980232e-08,2.980232e-08,
	2.980232e-08,2.980232e-08,0.004988182,2.980232e-08,0.00495712,2.980232e-08,0.007826605,0.002538008,2.980232e-08,0.002917703,0.003391779,2.980232e-08,0.009305687,2.980232e-08,2.980232e-08,1.891728e-08,9.624347e-05,0.0099455,0.006561823,
	0.008582851,2.980232e-08,0.008440943,0.006461001,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,2.980232e-08,0.005128652,0.009732047,2.980232e-08,0.002760089,2.980232e-08,2.980232e-08,2.980232e-08,0.00739862,
	0.007434984,2.980232e-08,0.0009944195,481.6896,0.00836362,2.980232e-08,0.008689271,2.980232e-08,0.003735749,2.980232e-08,0.00951215,0.003202765,2.980232e-08)

genout_m1=GenMatch(fit.func=lower.bound, Tr=college, X=Xmat, BalanceMatrix=out, estimand='ATT', M = 3, pop.size=1, max.generations=1, wait.generations=1, hard.generation.limit=F, starting.values=starts, nboots=0, ties=T, MemoryMatrix=F)	

1973

mout_m1=Match(Y=yppnscal,Tr=college, X=Xmat, estimand='ATT', M = 3, ties=T, Weight.matrix=genout_m1) 
summary(mout_m1)

matbal_m1 = MatchBalance(college ~ yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_m1,nboots=1000)
matbal_all1=MatchBalance(college ~ out, match.out = mout_m1)
mb1=percent.bal(matbal_all1)
mb1[[1]] # 0.4084507 :: 0.7098592

# Sensitivity Test

gamma_m1=hl.rbound(mout_m1, gamma=c(2,3,4,5), pr=.01, Y=yppnscal, paired=TRUE)

1982

mout_m4=Match(Y=y1982yppnscal[!is.na(y1982yppnscal)],Tr=college[!is.na(y1982yppnscal)], X=Xmat[!is.na(y1982yppnscal),], estimand='ATT', M = 3, ties=T, Weight.matrix=genout_m1) 
summary(mout_m4)

matbal_m4 = MatchBalance(college ~ yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_m4,nboots=1000)
matbal_all4=MatchBalance(college ~ out, match.out = mout_m4)
mb4=percent.bal(matbal_all4)
mb4[[1]] # 0.4084507 :: 0.5887324

# Sensitivity Test

gamma_m4=hl.rbound(mout_m4, gamma=c(2,3,4,5), pr=.01, Y=y1982yppnscal[!is.na(y1982yppnscal)], paired=TRUE)



#ATC: lin 1 to 3

matTemp=as.matrix(original)%*%diag(indicator[3,])   #469561
temp=list()
for(i in 1:81){
	if(indicator[3,i]==1){
		temp[[i]]=factor(matTemp[,i])
	}
	else{
		temp[[i]]=matTemp[,i]
	}
}

pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+
	temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+
	temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+
	temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
etahat1=pscore$linear.predictors


matTemp=as.matrix(original)%*%diag(indicator[4,])    #585646
temp=list()
for(i in 1:81){
	if(indicator[4,i]==1){
		temp[[i]]=factor(matTemp[,i])
	}
	else{
		temp[[i]]=matTemp[,i]
	}
}

pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+
	temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+
	temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+
	temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
etahat2=pscore$linear.predictors

Xmat=NULL
Xmat=cbind(etahat1,out)
for(i in 2:ncol(Xmat)){
	Xmat[,i]=lm(Xmat[,i]~etahat1)$residuals
}

Xmat=cbind(etahat2,Xmat)
for(i in 3:ncol(Xmat)){
	Xmat[,i]=lm(Xmat[,i]~etahat2)$residuals
}

starts=c(100,100.0071,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,739.1417,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,1.490116e-08,0.004539034,1.490116e-08,1.490116e-08,0.01517012,0.02992871,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.007834871,1.490116e-08,
	1.490116e-08,0.008113991,0.005919858,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.003814568,0.02419358,1.490116e-08,0.004334319,1.490116e-08,0.01948788,1.490116e-08,1.490116e-08,44.08746,
	1.490116e-08,1.490116e-08,13.53873,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.02696651,1.490116e-08,1.490116e-08,0.006232587,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.004892858,
	1.490116e-08,0.01903991,1.490116e-08,0.02379703,1.490116e-08,1.490116e-08,1.490116e-08,0.0034652,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,0.01299434,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.008628191,1.490116e-08,0.01575408,1.490116e-08,1.490116e-08,
	1.490116e-08,0.02771047,0.008514462,1.490116e-08,1.490116e-08,1.490116e-08,0.01750806,1.490116e-08,0.01683012,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.02491207,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,354.155,
	1.490116e-08,0.001847714,0.002439044,1.490116e-08,1.490116e-08,1.490116e-08,0.006686551,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,310.5782,1.490116e-08,1.490116e-08,0.01970852,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,0.007031653,1.490116e-08,1.490116e-08,0.01494422,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.007390262,1.490116e-08,0.006955788,1.490116e-08,0.0103705,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.0001336209,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.005842602,1.490116e-08,0.03062162,1.490116e-08,1.490116e-08,0.02512913,
	1.490116e-08,0.005845255,0.02964625,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.01362725,0.00318443,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.02585323,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.006876822,39.75011,0.0495586,1.490116e-08,1.490116e-08,0.01761491,0.02573231,
	1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.01040253,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.007306871,
	0.0003473068,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.009196998,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.01817766,
	0.001553622,3.720721e-06,89.98467,462.9543,210.3743,1.490116e-08,1.490116e-08,0.01362846,1.490116e-08,159.0117,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	0.008644248,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.008629846,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,0.001467478,573.9009,1.490116e-08,1.490116e-08,0.0001304165,1.490116e-08,1.490116e-08,
	1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,
	1.490116e-08,0.005418747,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08,1.490116e-08)

genout_m2=GenMatch(fit.func=lower.bound, Tr=college, X=Xmat, BalanceMatrix=out, estimand='ATC', M = 3, pop.size=1, max.generations=1, wait.generations=1, hard.generation.limit=F, starting.values=starts, nboots=0, ties=T, MemoryMatrix=F)	

1973

mout_m2=Match(Y=yppnscal,Tr=college, X=Xmat, estimand='ATC', M = 3, ties=T, Weight.matrix=genout_m2) 
summary(mout_m2)

matbal_m2 = MatchBalance(college ~ yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_m2,nboots=1000)
matbal_all2=MatchBalance(college ~ out, match.out = mout_m2)
mb2=percent.bal(matbal_all2)
mb2[[1]] # 0.4084507 :: 0.9014085

# Sensitivity Test

gamma_m2=hl.rbound(mout_m2, gamma=c(2,3,4,5), pr=.01, Y=yppnscal, paired=TRUE)

1982

mout_m5=Match(Y=y1982yppnscal[!is.na(y1982yppnscal)],Tr=college[!is.na(y1982yppnscal)], X=Xmat[!is.na(y1982yppnscal),], estimand='ATC', M = 3, ties=T, Weight.matrix=genout_m2) 
summary(mout_m5)

matbal_m5 = MatchBalance(college ~ yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_m5,nboots=1000)
matbal_all5=MatchBalance(college ~ out, match.out = mout_m5)
mb5=percent.bal(matbal_all5)
mb5[[1]] # 0.4084507 :: 0.7380282

# Sensitivity Test

gamma_m5=hl.rbound(mout_m5, gamma=c(2,3,4,5), pr=.01, Y=y1982yppnscal[!is.na(y1982yppnscal)], paired=TRUE)
  

#END