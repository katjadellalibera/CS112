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
# This file takes 'data_all_waves.dta' and creates 'data_rep', 'original', 'factors', 'out', and 'dat' data objects.  These
# objects are the global objects that the various R files require in order to compute the tables, statistics, and figures
# for the paper.  
#
# The 'data_all_waves.dta' is essentially the same data matrix as that provided by Kam and Palmer (2008), with additional 
# covariates included for the purposes of facilitating genetic matching analysis.  For replication purposes, we follow
# Kam and Palmer (2008) and thus also factorize the covariates for all of the propensity score analyses in 'Who Matches'.
  

# 2. Data

library(foreign)

data_rep=read.dta("data_all_waves.dta")
attach(data_rep)

# Factorizing Covariates

data_rep$yBlack=rep(0,length(yRace))
data_rep$yBlack[yRace==2]=1
data_rep$yWhite=rep(0,length(yRace))
data_rep$yWhite[yRace==1]=1
data_rep$pRep=pPID
data_rep$pRep[pPID==1]=0
data_rep$pRep[pPID==5|pPID==6|pPID==7]=1
data_rep$pRep[data_rep$pRep!=1]=0
data_rep$yRep=yPID
data_rep$yRep[yPID==1]=0
data_rep$yRep[yPID==5|yPID==6|yPID==7]=1
data_rep$yRep[data_rep$yRep!=1]=0
data_rep$pHHCollege=pEducHH
data_rep$pHHCollege[pEducHH>=5]=1
data_rep$pHHCollege[pEducHH<=4]=0
data_rep$pWCollege=pEducW
data_rep$pWCollege[pEducW>=5]=1
data_rep$pWCollege[pEducW<=4]=0
data_rep$pMednInc=pHHInc
data_rep$pMednInc[pHHInc>=8]=1
data_rep$pMednInc[pHHInc<=7]=0
data_rep$y1982yppnscal=(y1982vote80+y1982meeting+y1982other+y1982button+y1982money+y1982communicate+y1982demonstrate+y1982community)
data_rep$y1982yppnscal[which(data_rep$y1982yppnscal>80)]=NA

yPubAffZ = factor(yPubAff)
yNewspaperZ = factor(yNewspaper)
yRadioZ = factor(yRadio)
yTVZ = factor(yTV)
yMagazineZ = factor(yMagazine)
yFamTalkZ = factor(yFamTalk)
yFrTalkZ = factor(yFrTalk)
yAdultTalkZ = factor(yAdultTalk)
yPIDZ = factor(yPID)
ySPIDZ = factor(ySPID)
yGovtOpinionZ = factor(yGovtOpinion)
yGovtCrookZ = factor(yGovtCrook)
yGovtWasteZ = factor(yGovtWaste)
yTrGovtZ = factor(yTrGovt)
yGovtSmartZ = factor(yGovtSmart)
yGovt4AllZ = factor(yGovt4All)
yLifeWishZ = factor(yLifeWish)
yGLuckZ = factor(yGLuck)
yFPlansZ = factor(yFPlans)
yWinArgZ = factor(yWinArg)
yStrOpinionZ = factor(yStrOpinion)
yMChangeZ = factor(yMChange)
yTrOthersZ = factor(yTrOthers)
yOthHelpZ = factor(yOthHelp)
yOthFairZ = factor(yOthFair)
ySenateZ = factor(ySenate)
yTitoZ = factor(yTito)
yCourtZ = factor(yCourt)
yGovernZ = factor(yGovern)
yCCampZ = factor(yCCamp)
yFDRZ = factor(yFDR)
yKnowledgeZ = factor(yKnowledge)
yNextSchZ = factor(yNextSch)
yGPAZ = factor(yGPA)
ySchOfficerZ = factor(ySchOfficer)
ySchPublishZ = factor(ySchPublish)
yHobbyZ = factor(yHobby)
ySchClubZ = factor(ySchClub)
yOccClubZ = factor(yOccClub)
yNeighClubZ = factor(yNeighClub)
yRelClubZ = factor(yRelClub)
yYouthOrgZ = factor(yYouthOrg)
yMiscClubZ = factor(yMiscClub)
yClubLevZ = factor(yClubLev)
yPhoneZ = factor(yPhone)
yGenZ = factor(yGen)
yRaceZ = factor(yRace)
yEgoAZ = factor(yEgoA)
yEgoBZ = factor(yEgoB)
yTrustZ = factor(yTrust)
yCynicZ = factor(yCynic)
pNewspaperZ = factor(pNewspaper)
pRadioZ = factor(pRadio)
pTVZ = factor(pTV)
pMagazineZ = factor(pMagazine)
pGLuckZ = factor(pGLuck)
pFPlansZ = factor(pFPlans)
pWinArgZ = factor(pWinArg)
pStrOpinionZ = factor(pStrOpinion)
pMChangeZ = factor(pMChange)
pTrOthersZ = factor(pTrOthers)
pOthHelpZ = factor(pOthHelp)
pOthFairZ = factor(pOthFair)
pPIDZ = factor(pPID)
pSPIDZ = factor(pSPID)
pVoteZ = factor(pVote)
pPersuadeZ = factor(pPersuade)
pRallyZ = factor(pRally)
pOthActZ = factor(pOthAct)
pPolClubZ = factor(pPolClub)
pButtonZ = factor(pButton)
pMoneyZ = factor(pMoney)
pActFrqZ = factor(pActFrq)
pGovtOpinionZ = factor(pGovtOpinion)
pGovtCrookZ = factor(pGovtCrook)
pGovtWasteZ = factor(pGovtWaste)
pTrGovtZ = factor(pTrGovt)
pGovtSmartZ = factor(pGovtSmart)
pGovt4AllZ = factor(pGovt4All)
pLifeWishZ = factor(pLifeWish)
pEmployZ = factor(pEmploy)
pEducHHZ = factor(pEducHH)
pEducWZ = factor(pEducW)
pChurchOrgZ = factor(pChurchOrg)
pFratOrgZ = factor(pFratOrg)
pProOrgZ = factor(pProOrg)
pCivicOrgZ = factor(pCivicOrg)
pCLOrgZ = factor(pCLOrg)
pNeighClubZ = factor(pNeighClub)
pSportClubZ = factor(pSportClub)
pInfClubZ = factor(pInfClub)
pFarmGrZ = factor(pFarmGr)
pWomenClubZ = factor(pWomenClub)
pMiscClubZ = factor(pMiscClub)
pClubLevZ = factor(pClubLev)
pFIncZ = factor(pFInc)
pHHIncZ = factor(pHHInc)
pOwnHomeZ = factor(pOwnHome)
pSenateZ = factor(pSenate)
pTitoZ = factor(pTito)
pCourtZ = factor(pCourt)
pGovernZ = factor(pGovern)
pCCampZ = factor(pCCamp)
pFDRZ = factor(pFDR)
pKnowledgeZ = factor(pKnowledge)
pGenZ = factor(pGen)
pRaceZ = factor(pRace)
pParticipate1Z=factor(pParticipate1)
pParticipate2Z=factor(pParticipate2)
pHHCollegeZ=factor(data_rep$pHHCollege)
pWCollegeZ=factor(data_rep$pWCollege)
pMednIncZ=factor(data_rep$pMednInc)


factors = data.frame(yPubAffZ, yNewspaperZ, yRadioZ, yTVZ, yMagazineZ, yFamTalkZ, yFrTalkZ, yAdultTalkZ, yPIDZ, ySPIDZ, yGovtOpinionZ, yGovtCrookZ, yGovtWasteZ, yTrGovtZ, yGovtSmartZ, yGovt4AllZ, yLifeWishZ, yGLuckZ, yFPlansZ, yWinArgZ, yStrOpinionZ, 
	yMChangeZ, yTrOthersZ, yOthHelpZ, yOthFairZ, ySenateZ, yTitoZ, yCourtZ, yGovernZ, yCCampZ, yFDRZ, yKnowledgeZ, yNextSchZ, yGPAZ, ySchOfficerZ, ySchPublishZ, yHobbyZ, ySchClubZ, yOccClubZ, yNeighClubZ, yRelClubZ, yYouthOrgZ, yMiscClubZ, yClubLevZ, 
	yPhoneZ, yGenZ, yRaceZ, yEgoAZ, yEgoBZ, yTrustZ, yCynicZ, pNewspaperZ, pRadioZ, pTVZ, pMagazineZ, pLifeWishZ, pGLuckZ, pFPlansZ, pWinArgZ, pStrOpinionZ, pMChangeZ, pTrOthersZ, pOthHelpZ, pOthFairZ, pPIDZ, pSPIDZ, pVoteZ, pPersuadeZ, pRallyZ, pOthActZ, 
	pPolClubZ, pButtonZ, pMoneyZ, pActFrqZ, pGovtOpinionZ, pGovtCrookZ, pGovtWasteZ, pTrGovtZ, pGovtSmartZ, pGovt4AllZ, pEmployZ, pEducHHZ, pEducWZ, pChurchOrgZ, pFratOrgZ, pProOrgZ, pCivicOrgZ, pCLOrgZ, pNeighClubZ, pSportClubZ, pInfClubZ, pFarmGrZ, 
	pWomenClubZ, pMiscClubZ, pClubLevZ, pFIncZ, pHHIncZ, pOwnHomeZ, pSenateZ, pTitoZ, pCourtZ, pGovernZ, pCCampZ, pFDRZ, pKnowledgeZ, pGenZ, pRaceZ, pParticipate1Z, pParticipate2Z,pHHCollegeZ,pWCollegeZ,pMednIncZ)
original = data.frame(cbind(yPubAffZ,yNewspaperZ , yRadioZ , yMagazineZ , yFamTalkZ , yFrTalkZ , yAdultTalkZ , ySPIDZ , yGovtOpinionZ , yGovtCrookZ , yGovtWasteZ , yTrGovtZ , yGovtSmartZ , yGovt4AllZ , yLifeWishZ , yGLuckZ , yFPlansZ , yWinArgZ , 
	yStrOpinionZ , yMChangeZ , yTrOthersZ , yOthHelpZ , yOthFairZ , yKnowledgeZ , yNextSchZ , yGPAZ , ySchOfficerZ , ySchPublishZ , yHobbyZ , ySchClubZ , yOccClubZ , yNeighClubZ , yRelClubZ , yYouthOrgZ , yClubLevZ , yPhoneZ , yGenZ , yRaceZ , pNewspaperZ , 
	pRadioZ , pTVZ , pMagazineZ , pLifeWishZ , pGLuckZ , pFPlansZ , pWinArgZ , pStrOpinionZ , pMChangeZ , pTrOthersZ , pOthHelpZ , pOthFairZ , pSPIDZ , pVoteZ , pPersuadeZ , pRallyZ , pOthActZ , pPolClubZ , pButtonZ , pMoneyZ , pGovtOpinionZ , pGovtCrookZ , 
	pGovtWasteZ , pTrGovtZ , pGovtSmartZ , pGovt4AllZ , pEmployZ , pEducHHZ , pChurchOrgZ , pFratOrgZ , pProOrgZ , pCivicOrgZ , pCLOrgZ , pNeighClubZ , pSportClubZ , pInfClubZ , pFarmGrZ , pWomenClubZ , pClubLevZ , pHHIncZ , pOwnHomeZ , pKnowledgeZ))


dat=data_rep[,13:121]
dat$pParticipate2=pParticipate2*5
dat$pParticipate1=pParticipate1*6
dat$yKnowledge=yKnowledge*6
dat$pKnowledge=pKnowledge*6
dat[dat<1.1 & dat > 0.9]=1
dat[dat<2.1 & dat > 1.9]=2
dat[dat<3.1 & dat > 2.9]=3
dat[dat<4.1 & dat > 3.9]=4
dat[dat<5.1 & dat > 4.9]=5
dat[dat<6.1 & dat > 5.9]=6
dat[dat<7.1 & dat > 6.9]=7

dummy=function(X=dat){
	vars=list()
	
	for(i in 1:ncol(X)){
		
		if(min(X[,i]==0)){
			X[,i]=X[,i] + 1
		}
		
		if(max(as.numeric(levels(factor(X[,i]))))==1){
			vars[[i]]=(X[,i])
		}
		
		else{
			temps=matrix(0,length(X[,i]),(max(as.numeric(levels(factor(X[,i]))))))  
			for(j in 1:(max(as.numeric(levels(factor(X[,i]))))))  {  
				temps[,j][X[,i]==j]=1
			}
			vars[[i]]=temps
		}
	}
	
	bound=function(X=vars){
		out=X[[1]]
		for(i in 2:length(X)){
			out=cbind(out,X[[i]])
		}
		return(out)
	}
	
	out=bound(vars)
}

out=NULL
out=dummy(dat)


# 3. Save

#save(data_rep,original,factors,out,dat,indicator,file="WhoMatches.Rdata")



# END