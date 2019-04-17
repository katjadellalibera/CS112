#######################
#   John Henderson    #            
#   JoP Who Matches   #
#                     #
#   Iteration Funcs   #
#                     #
#   Sept. 12, 2009    #
#######################



# PSCORE ITERATION FUNCTION

# Function to Implement Random Permuted Propensity Score Estimation Procedure
#
# Author: John Henderson
# Date: 7/12/2008
#
#
# Takes an indicator matrix (matrix of indicator vectors) and accordingly permutes from a matrix of covariates, 
#  to estimate a propensity score model m number of times. Ultimately estimates 6 models:
#
#  (a) ATT: 1 to 1, predicted probability
#  (b) ATT: 1 to 3, predicted probability
#  (c) ATT: 1 to 1, linear predictor
#  (d) ATT: 1 to 3, linear predictor
#  (e) ATC: 1 to 3, linear predictor
#  (f) ATC: 1 to 1, linear predictor

# indicator - A matrix of 1 and 0 with number of rows, m, equaling the number of iterations desired, and the number of 
#             columns the number of covariates being sampled 
   

pscore.eval=function(indicator){
	
	out=out
	indicator=as.numeric(indicator)
	original=original
	data_rep=data_rep
	yppnscal=data_rep$yppnscal
	college=data_rep$college
	
	
	matTemp=as.matrix(original)%*%diag(indicator)
	temp=list()
	for(i in 1:81){
		if(indicator[i]==1){
			temp[[i]]=factor(matTemp[,i])
		}
		else{
			temp[[i]]=matTemp[,i]
		}
	}
	
	pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
	
	#Predicted Probabilities 
	#ATT: 1 to 1	
	
	etahat=pscore$fitted.values
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATT", M = 1,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects1=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))
	
	#Predicted Probabilities 
	#ATT: 1 to 3	
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATT", M = 3,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects2=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))
	
	#Linear Predictors  
	#ATT: 1 to 1		
	
	etahat=pscore$linear.predictors
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATT", M = 1,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects3=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))
	
	#Linear Predictors 
	#ATT: 1 to 3	
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATT", M = 3,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects4=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))
	
	#Linear Predictors 
	#ATC: 1 to 1
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATC", M = 1,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects5=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))

	#Linear Predictors 
	#ATC: 1 to 3	
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand="ATC", M = 3,ties=TRUE) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
#	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object3=0
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yRace1[mat$index.treated],data_rep$yRace1[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRace2[mat$index.treated],data_rep$yRace2[mat$index.control],mat$weights)
	object8=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
#	object9=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=0
	object10=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object11=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object12=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object13=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
#	object14=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=0
	object15=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object16=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
#	object17=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=0
	object18=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
#	object19=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=0
	object20=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
#	object21=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=0
	object22=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object23=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object24=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
#	object25=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=0
	object26=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object27=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/24)
	object28=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/24
	object29=var(objects.a)
	
	object30=(length(which(outs>.1)))/length(outs) 
	object31=mean(outs) 
	object32=var(outs)
	
	object33=(length(which(etahat>.99))/1254)
	object34=(length(which(etahat>.95))/1254)
	object35=(length(which(etahat>.9))/1254)
	object36=(length(which(etahat>.75))/1254)
	object37=(length(which(etahat>.5))/1254)
	object38=(length(which(etahat>.25))/1254)
	object39=(length(which(etahat>.1))/1254)
	
	objects6=cbind(c(object1,object2,objects.a,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39))
	
	objects=rbind(objects1,objects2,objects3,objects4,objects5,objects6)
	
	return(objects)
}



# GENMATCH ITERATION FUNCTION

# Function to Implement Lower Bounded Genetic Matching Algorithm Starting at the Best Propensity Score Estimate 
#
# Author: John Henderson
# Date: 5/12/2009
#
#
# Takes an indicator vector associated with the best propensity score model, and using that model as the starting point 
#  in a GenMatch run, that either maximizes setting lower bounds or maximizes the average proportion of balanced variables
#  at the p > .1 level
#
# indicator - A vector of 1 and 0
# estimand - Either ATC or ATT
# M - The number of control observations to be matched to each treated observation
# pop.size, max.gen, wait.gen - GenMatch quantities specifying the number of multivariate matching weights to be used to 
#              optimize the fit function
# linear - Specifies either 'linpred' or 'predprob', indicating whether the linear predictor or predicted probability will
#          be used for the propensity score starting point
#
# See http://sekhon.berkeley.edu/matching/GenMatch.html for more details. 


gen.eval=function(indicator, estimand, M, pop.size, max.gen, wait.gen,linear='linpred'){
	
	out=out
	indicator=as.numeric(indicator)
	original=original
	data_rep=data_rep
	yppnscal=data_rep$yppnscal
	college=data_rep$college

	
	matTemp=as.matrix(original)%*%diag(indicator)
	temp=list()
	for(i in 1:81){
		if(indicator[i]==1){
			temp[[i]]=factor(matTemp[,i])
		}
		else{
			temp[[i]]=matTemp[,i]
		}
	}
	
	pscore=glm(college~temp[[1]]+temp[[2]]+temp[[3]]+temp[[4]]+temp[[5]]+temp[[6]]+temp[[7]]+temp[[8]]+temp[[9]]+temp[[10]]+temp[[11]]+temp[[12]]+temp[[13]]+temp[[14]]+temp[[15]]+temp[[16]]+temp[[17]]+temp[[18]]+temp[[19]]+temp[[20]]+temp[[21]]+temp[[22]]+temp[[23]]+temp[[24]]+temp[[25]]+temp[[26]]+temp[[27]]+temp[[28]]+temp[[29]]+temp[[30]]+temp[[31]]+temp[[32]]+temp[[33]]+temp[[34]]+temp[[35]]+temp[[36]]+temp[[37]]+temp[[38]]+temp[[39]]+temp[[40]]+temp[[41]]+temp[[42]]+temp[[43]]+temp[[44]]+temp[[45]]+temp[[46]]+temp[[47]]+temp[[48]]+temp[[49]]+temp[[50]]+temp[[51]]+temp[[52]]+temp[[53]]+temp[[54]]+temp[[55]]+temp[[56]]+temp[[57]]+temp[[58]]+temp[[59]]+temp[[60]]+temp[[61]]+temp[[62]]+temp[[63]]+temp[[64]]+temp[[65]]+temp[[66]]+temp[[67]]+temp[[68]]+temp[[69]]+temp[[70]]+temp[[71]]+temp[[72]]+temp[[73]]+temp[[74]]+temp[[75]]+temp[[76]]+temp[[77]]+temp[[78]]+temp[[79]]+temp[[80]]+temp[[81]],family=binomial(link=logit))  
	
	if(linear=='linpred'){
		etahat=pscore$linear.predictors
	}
	else if(linear=='predprob'){
		etahat=pscore$fitted.values
	}
	
	#Pscore Match
	
	mat=Match(Y=yppnscal,Tr=college, X=etahat, estimand=estimand, M = M, ties=T) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yBlack[mat$index.treated],data_rep$yBlack[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
	object8=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object10=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object11=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object12=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
	object13=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object15=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
	object16=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
	object18=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
	object20=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object22=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object23=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
	object24=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	object26=ksboot(data_rep$yMagazine[mat$index.treated],data_rep$yMagazine[mat$index.control])$ks.boot.pvalue
	object27=ttest(data_rep$yMagazine[mat$index.treated],data_rep$yMagazine[mat$index.control],mat$weights)
	object28=ksboot(data_rep$ySchClub[mat$index.treated],data_rep$ySchClub[mat$index.control])$ks.boot.pvalue
	object29=ttest(data_rep$ySchClub[mat$index.treated],data_rep$ySchClub[mat$index.control],mat$weights)
	object30=ksboot(data_rep$yFrTalk[mat$index.treated],data_rep$yFrTalk[mat$index.control])$ks.boot.pvalue
	object31=ttest(data_rep$yFrTalk[mat$index.treated],data_rep$yFrTalk[mat$index.control],mat$weights)
	object32=ksboot(data_rep$yNeighClub[mat$index.treated],data_rep$yNeighClub[mat$index.control])$ks.boot.pvalue
	object33=ttest(data_rep$yNeighClub[mat$index.treated],data_rep$yNeighClub[mat$index.control],mat$weights)
	object34=ksboot(data_rep$pTrOthers[mat$index.treated],data_rep$pTrOthers[mat$index.control])$ks.boot.pvalue
	object35=ttest(data_rep$pTrOthers[mat$index.treated],data_rep$pTrOthers[mat$index.control],mat$weights)
	object36=ksboot(data_rep$pGovtWaste[mat$index.treated],data_rep$pGovtWaste[mat$index.control])$ks.boot.pvalue
	object37=ttest(data_rep$pGovtWaste[mat$index.treated],data_rep$pGovtWaste[mat$index.control],mat$weights)
	object38=ksboot(data_rep$pChurchOrg[mat$index.treated],data_rep$pChurchOrg[mat$index.control])$ks.boot.pvalue
	object39=ttest(data_rep$pChurchOrg[mat$index.treated],data_rep$pChurchOrg[mat$index.control],mat$weights)
	object40=ksboot(data_rep$pProOrg[mat$index.treated],data_rep$pProOrg[mat$index.control])$ks.boot.pvalue
	object41=ttest(data_rep$pProOrg[mat$index.treated],data_rep$pProOrg[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39,object40,object41)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs=matrix(c(apply(out,2,bal.vec)))
	
	object42=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/40)
	object43=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/40
	object44=var(objects.a)
	
	object45=(length(which(outs>.1)))/length(outs) 
	object46=mean(outs) 
	object47=var(outs)
	
	object48=(length(which(etahat>.99))/1254)
	object49=(length(which(etahat>.95))/1254)
	object50=(length(which(etahat>.9))/1254)
	object51=(length(which(etahat>.75))/1254)
	object52=(length(which(etahat>.5))/1254)
	object53=(length(which(etahat>.25))/1254)
	object54=(length(which(etahat>.1))/1254)
	
	objects1=cbind(c(object1,object2,objects.a,object42,object43,object44,object45,object46,object47,object48,object49,object50,object51,object52,object53,object54))
	
	Xmat=cbind(etahat,out)
	for(i in 2:ncol(Xmat)){
		Xmat[,i]=lm(Xmat[,i]~etahat)$residuals
	}
	
	#Xmat=cbind(etahat,out[,which(outs<.1)])
	#for(i in 2:ncol(Xmat)){
	#	Xmat[,i]=lm(Xmat[,i]~etahat)$residuals
	#}
	
	starts=rep(0,ncol(Xmat))
	starts[1]=100

	#GenMatch
	
	genout=GenMatch(fit.func=lower.bound, Tr=college, X=Xmat, BalanceMatrix=out, estimand=estimand, M = M, pop.size=pop.size, max.generations=max.gen, wait.generations=wait.gen, hard.generation.limit=F, starting.values=starts, nboots=0, ties=T, MemoryMatrix=F)	
	
	mat=Match(Y=yppnscal,Tr=college, X=Xmat, estimand=estimand, M = M, ties=T, Weight.matrix=genout) 
	
	object1=mat$est	
	object2=((1-pnorm(abs(mat$est/mat$se)))*2)
	object3=ksboot(data_rep$yGPA[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object4=ttest(data_rep$yGPA[mat$index.treated],data_rep$yGPA[mat$index.control],mat$weights)
	object5=ttest(data_rep$yGen[mat$index.treated],data_rep$yGen[mat$index.control],mat$weights)
	object6=ttest(data_rep$yBlack[mat$index.treated],data_rep$yBlack[mat$index.control],mat$weights)
	object7=ttest(data_rep$yRep[mat$index.treated],data_rep$yRep[mat$index.control],mat$weights)
	object8=ksboot(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control])$ks.boot.pvalue
	object9=ttest(data_rep$yKnowledge[mat$index.treated],data_rep$yKnowledge[mat$index.control],mat$weights)
	object10=ttest(data_rep$yNextSch[mat$index.treated],data_rep$yNextSch[mat$index.control],mat$weights)
	object11=ttest(data_rep$pVote[mat$index.treated],data_rep$pVote[mat$index.control],mat$weights)
	object12=ttest(data_rep$pPersuade[mat$index.treated],data_rep$pPersuade[mat$index.control],mat$weights)
	object13=ksboot(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control])$ks.boot.pvalue
	object14=ttest(data_rep$pParticipate2[mat$index.treated],data_rep$pParticipate2[mat$index.control],mat$weights)
	object15=ttest(data_rep$pEmploy[mat$index.treated],data_rep$pEmploy[mat$index.control],mat$weights)
	object16=ksboot(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control])$ks.boot.pvalue
	object17=ttest(data_rep$pEducHH[mat$index.treated],data_rep$pEducHH[mat$index.control],mat$weights)
	object18=ksboot(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control])$ks.boot.pvalue
	object19=ttest(data_rep$pEducW[mat$index.treated],data_rep$pEducW[mat$index.control],mat$weights)
	object20=ksboot(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control])$ks.boot.pvalue
	object21=ttest(data_rep$pHHInc[mat$index.treated],data_rep$pHHInc[mat$index.control],mat$weights)
	object22=ttest(data_rep$pOwnHome[mat$index.treated],data_rep$pOwnHome[mat$index.control],mat$weights)
	object23=ttest(data_rep$pRep[mat$index.treated],data_rep$pRep[mat$index.control],mat$weights)
	object24=ksboot(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control])$ks.boot.pvalue
	object25=ttest(data_rep$pKnowledge[mat$index.treated],data_rep$pKnowledge[mat$index.control],mat$weights)
	object26=ksboot(data_rep$yMagazine[mat$index.treated],data_rep$yMagazine[mat$index.control])$ks.boot.pvalue
	object27=ttest(data_rep$yMagazine[mat$index.treated],data_rep$yMagazine[mat$index.control],mat$weights)
	object28=ksboot(data_rep$ySchClub[mat$index.treated],data_rep$ySchClub[mat$index.control])$ks.boot.pvalue
	object29=ttest(data_rep$ySchClub[mat$index.treated],data_rep$ySchClub[mat$index.control],mat$weights)
	object30=ksboot(data_rep$yFrTalk[mat$index.treated],data_rep$yFrTalk[mat$index.control])$ks.boot.pvalue
	object31=ttest(data_rep$yFrTalk[mat$index.treated],data_rep$yFrTalk[mat$index.control],mat$weights)
	object32=ksboot(data_rep$yNeighClub[mat$index.treated],data_rep$yNeighClub[mat$index.control])$ks.boot.pvalue
	object33=ttest(data_rep$yNeighClub[mat$index.treated],data_rep$yNeighClub[mat$index.control],mat$weights)
	object34=ksboot(data_rep$pTrOthers[mat$index.treated],data_rep$pTrOthers[mat$index.control])$ks.boot.pvalue
	object35=ttest(data_rep$pTrOthers[mat$index.treated],data_rep$pTrOthers[mat$index.control],mat$weights)
	object36=ksboot(data_rep$pGovtWaste[mat$index.treated],data_rep$pGovtWaste[mat$index.control])$ks.boot.pvalue
	object37=ttest(data_rep$pGovtWaste[mat$index.treated],data_rep$pGovtWaste[mat$index.control],mat$weights)
	object38=ksboot(data_rep$pChurchOrg[mat$index.treated],data_rep$pChurchOrg[mat$index.control])$ks.boot.pvalue
	object39=ttest(data_rep$pChurchOrg[mat$index.treated],data_rep$pChurchOrg[mat$index.control],mat$weights)
	object40=ksboot(data_rep$pProOrg[mat$index.treated],data_rep$pProOrg[mat$index.control])$ks.boot.pvalue
	object41=ttest(data_rep$pProOrg[mat$index.treated],data_rep$pProOrg[mat$index.control],mat$weights)
	
	objects.a=rbind(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26,object27,object28,object29,object30,object31,object32,object33,object34,object35,object36,object37,object38,object39,object40,object41)
	
	bal.vec=function(out){
		mat=mat
		out=out
		outs=matrix()
		outs=ttest(out[mat$index.treated],out[mat$index.control],mat$weights)
		return(outs)
	}
	
	outs1=matrix(c(apply(out,2,bal.vec)))
	
	object42=(length(which(c(object3,object4,object5,object6,object7,object8,object9,object10,object11,object12,object13,object14,object15,object16,object17,object18,object19,object20,object21,object22,object23,object24,object25,object26)>.1))/40)
	object43=(object3+object4+object5+object6+object7+object8+object9+object10+object11+object12+object13+object14+object15+object16+object17+object18+object19+object20+object21+object22+object23+object24+object25+object26)/40
	object44=var(objects.a)
	
	object45=(length(which(outs1>.1)))/length(outs1) 
	object46=mean(outs1) 
	object47=var(outs1)
	
	object48=(length(which(etahat>.99))/1254)
	object49=(length(which(etahat>.95))/1254)
	object50=(length(which(etahat>.9))/1254)
	object51=(length(which(etahat>.75))/1254)
	object52=(length(which(etahat>.5))/1254)
	object53=(length(which(etahat>.25))/1254)
	object54=(length(which(etahat>.1))/1254)
	
	objects2=cbind(c(object1,object2,objects.a,object42,object43,object44,object45,object46,object47,object48,object49,object50,object51,object52,object53,object54))
	
	objects=cbind(objects1,objects2)
	out.objects=list()
	out.objects[[1]]=objects
	out.objects[[2]]=genout
	return(out.objects)
}



# END