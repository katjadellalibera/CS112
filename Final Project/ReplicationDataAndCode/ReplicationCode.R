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
# This is the code for the replication of Kam and Palmers (2008) propensity score matching analysis. This code is 
# slightly modified and expanded from the original code provided by Kam and Palmer (2008), and correspondes to their article 
# 'Reconsidering the Effects of Education on Political Participation'.  
# 
# In addition to the authors propensity score, included below are MatchBalance statistics, placbeo tests, and alternative 
# matching models.
#  
# Tables 1 - 4 in the paper "Who Matches?", are produced below.  These are the models 'mout_r1', 'mout_r1D', 'mout_r3', 
# and 'mout_r3D'.  Balance statistisc are found in 'matbal_r1' and 'matbal_r3'.
#
# Additionally, diagnostic statistics are also estimated below.  These includes the moments of the Kam and Palmer (2008) 
# propesity score, as well as its degree of overlap, among other things.
 

library(MASS)
library(Matching)
 
source('funcs.R')
source('objects.R')   
#load("WhoMatches.Rdata")
attach(data_rep)
attach(factors)

# KP Pscore
#  Note the the Z suffix indicates that the covariate is factorized following Kam and Palmer (2008)

model=formula(college ~ yPubAffZ + yNewspaperZ + yRadioZ + yMagazineZ + yFamTalkZ + yFrTalkZ + yAdultTalkZ + ySPIDZ + yGovtOpinionZ + yGovtCrookZ + yGovtWasteZ + yTrGovtZ + yGovtSmartZ + yGovt4AllZ + yLifeWishZ + yGLuckZ + yFPlansZ + yWinArgZ + yStrOpinionZ + 
	yMChangeZ + yTrOthersZ + yOthHelpZ + yOthFairZ + yKnowledgeZ + yNextSchZ + yGPAZ + ySchOfficerZ + ySchPublishZ + yHobbyZ + ySchClubZ + yOccClubZ + yNeighClubZ + yRelClubZ + yYouthOrgZ + yClubLevZ + yPhoneZ + yGenZ + yRaceZ + pNewspaperZ + pRadioZ + pTVZ + 
	pMagazineZ + pLifeWishZ + pGLuckZ + pFPlansZ + pWinArgZ + pStrOpinionZ + pMChangeZ + pTrOthersZ + pOthHelpZ + pOthFairZ + pSPIDZ + pVoteZ + pPersuadeZ + pRallyZ + pOthActZ + pPolClubZ + pButtonZ + pMoneyZ + pGovtOpinionZ + pGovtCrookZ + pGovtWasteZ + pTrGovtZ + 
	pGovtSmartZ + pGovt4AllZ + pEmployZ + pEducHHZ + pChurchOrgZ + pFratOrgZ + pProOrgZ + pCivicOrgZ + pCLOrgZ + pNeighClubZ + pSportClubZ + pInfClubZ + pFarmGrZ + pWomenClubZ + pClubLevZ + pHHIncZ + pOwnHomeZ + pKnowledgeZ)
pscore=glm(model,family=binomial(link=logit))
etahat=pscore$fitted.values


# Propensity Score Matching
#  Note that 'out' is a n x 355 matrix, which is produced by dichotomizing 109 (86) covariates.  

#1973: 1 to 3

mout_r1 = Match(Y = yppnscal, Tr = college, X = etahat, estimand="ATT", M = 3) 
summary(mout_r1)

matbal_r1 = MatchBalance(college ~ etahat + yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_r1,nboots=1000)
matbal_all1=MatchBalance(college ~ out, match.out = mout_r1)
mb1=percent.bal(matbal_all1)
mb1[[1]] # 0.4084507 :: 0.2535211

# Outliers Dropped: 723, 676, 1061, 337, 595

mout_r1D = Match(Y= yppnscal[-c(723,676,1061,337,595)], Tr = college[-c(723,676,1061,337,595)], X = etahat[-c(723,676,1061,337,595)], estimand="ATT", M = 3) 
summary(mout_r1D)


#1973: 1 to 1

mout_r2 = Match(Y=yppnscal, Tr=college, X=etahat, estimand="ATT", M=1)   
summary(mout_r2)

matbal_r2 = MatchBalance(college ~ etahat + yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_r2,nboots=1000)
matbal_all2=MatchBalance(college ~ out, match.out = mout_r2)
mb2=percent.bal(matbal_all2)
mb2[[1]] # 0.4084507 :: 0.2366197

# Outliers Dropped: 723, 676, 1061, 337, 595

mout_r2D = Match(Y=yppnscal[-c(723,676,1061,337,595)], Tr=college[-c(723,676,1061,337,595)], X=etahat[-c(723,676,1061,337,595)], estimand="ATT", M=1)   
summary(mout_r2D)


#1982: 1 to 3

mout_r3 = Match(Y = y1982yppnscal[!is.na(y1982yppnscal)], Tr = college[!is.na(y1982yppnscal)], X = etahat[!is.na(y1982yppnscal)], estimand="ATT", M = 3) 
summary(mout_r3)

matbal_r3 = MatchBalance(college ~ etahat + yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_r3,nboots=1000)
matbal_all3=MatchBalance(college ~ out, match.out = mout_r3)
mb3=percent.bal(matbal_all3)
mb3[[1]] # 0.4084507 :: 0.2394366

# Outliers Dropped [after dropping missings in 1982 and re-indexing]:  613, 572, 897, 337, 507 

mout_r3D = Match(Y = y1982yppnscal[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], Tr = college[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], X = etahat[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], estimand="ATT", M = 3) 
summary(mout_r3D)


#1982: 1 to 1

mout_r4 = Match(Y = y1982yppnscal[!is.na(y1982yppnscal)], Tr = college[!is.na(y1982yppnscal)], X = etahat[!is.na(y1982yppnscal)], estimand="ATT", M = 1) 
summary(mout_r4)

matbal_r4 = MatchBalance(college ~ etahat + yGPA + yGen + yBlack + yRep+ yKnowledge + yNextSch + pVote + pPersuade + pParticipate2 + pEmploy + pEducHH + pEducW + pHHInc + pOwnHome + pRep  + pKnowledge, match.out = mout_r4,nboots=1000)
matbal_all4=MatchBalance(college ~ out, match.out = mout_r4)
mb4=percent.bal(matbal_all4)
mb4[[1]] # 0.4084507 :: 0.2732394

# Outliers Dropped [after dropping missings in 1982 and re-indexing]:  613, 572, 897, 337, 507 

mout_r4D = Match(Y = y1982yppnscal[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], Tr = college[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], X = etahat[!is.na(y1982yppnscal)][-c(572,897,296,613,507)], estimand="ATT", M = 1) 
summary(mout_r4D)


# Pscore Stats/Diagnostics

# Control Outlers Index (1973): 723, 676, 1061, 337, 595
# Control Outlers Index (1982): 613, 572, 897, 337, 507 

# Moments 

caliper = .25
mean(etahat)                                                                # 0.6403509
mean(etahat[college==1])                                                    # 0.8338374
mean(etahat[college==0])                                                    # 0.2958504
sd(etahat)                                                                  # 0.3512628
sd(etahat)*caliper                                                          # 0.0878157

# Max/Min/95th Percentile

sort(etahat[college==0])[round(.95*length(etahat[college==0]))+1]           # 0.8280431 
sort(etahat[college==1])[round(.95*length(etahat[college==1]))+1]           # 0.9998174 
max(etahat[college==1])                                                     # 0.9999998
min(etahat[college==1])                                                     # 0.02330369
max(etahat[college==0])                                                     # 0.9889442
min(etahat[college==0])                                                     # 1.196065e-10

# No Overlap at the Endpoints

length(which(etahat[college==0]>.9))/length(etahat[college==0])             # 0.03104213
length(which(etahat[college==0]>.95))/length(etahat[college==0])            # 0.01108647
length(which(etahat[college==1]>.9))/length(etahat[college==0])             # 0.5728518
length(which(etahat[college==1]>.95))/length(etahat[college==0])            # 0.4458281

length(which(etahat[college==1] > max(etahat[college==0])))/nrow(data_rep)  # 0.1714514
length(which(etahat[college==0] < min(etahat[college==1])))/nrow(data_rep)  # 0.05582137

# Outliers Index (1973)

which(etahat[college==0]>.95) 
sort(etahat[c(337,595,676,723,1061)])
index_r1=cbind(mout_r1$index.treated,mout_r1$index.control)

# Corresponding Matches of Outliers to College Attenders (Treatment)

length(which(index_r1[,2]==337))                                            # 60
length(which(index_r1[,2]==595))                                            # 58
length(which(index_r1[,2]==676))                                            # 316
length(which(index_r1[,2]==723))                                            # 287
length(which(index_r1[,2]==1061))                                           # 335

length(which(index_r1[,2]==337))/nrow(index_r1)                             # 0.02451982
length(which(index_r1[,2]==595))/nrow(index_r1)                             # 0.02370249
length(which(index_r1[,2]==676))/nrow(index_r1)                             # 0.1291377
length(which(index_r1[,2]==723))/nrow(index_r1)                             # 0.1172865
length(which(index_r1[,2]==1061))/nrow(index_r1)                            # 0.1369023

(length(which(index_r1[,2]==337))+length(which(index_r1[,2]==595))+length(which(index_r1[,2]==676))+length(which(index_r1[,2]==723))+length(which(index_r1[,2]==1061)))/nrow(index_r1) # 0.4315488

# Participatory Outliers: 1973 & 1982

mean(yppnscal[college==1 & etahat>.95])                                     # 3.26257
mean(yppnscal[college==0 & etahat>.95])                                     # 4
mean(yppnscal[college==1])                                                  # 2.793275
mean(yppnscal[college==0])                                                  # 1.427938

mean(data_rep$y1982yppnscal[college==1 & etahat>.95],na.rm=T)               # 3.319355
mean(data_rep$y1982yppnscal[college==0 & etahat>.95],na.rm=T)               # 4.4
mean(data_rep$y1982yppnscal[college==1],na.rm=T)                            # 3.041543
mean(data_rep$y1982yppnscal[college==0],na.rm=T)                            # 1.976064



# END