### 1 ###
## a ##
library(Matching)
data("lalonde")
attach(lalonde)

X<- cbind(age,educ,black, hisp, married, nodegr, u74, u75, re75, re74)
genout <- GenMatch(Tr=treat, X=X)
mout<- Match(Y=re78,Tr=treat,X=X, Weight.matrix = genout)
summary(mout)
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)

## b ## 
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=1,
                   pop.size=16, max.generations=10, wait.generations=1)
#The outcome variable
Y=re78/1000
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
mout <- Match(Y=Y, Tr=treat, X=X, Weight.matrix=genout)
summary(mout)
#Let's determine if balance has actually been obtained on the variables of interest
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)

## c ##
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat)
#The outcome variable
Y=re78/1000
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
mout <- Match(Y=Y, Tr=treat, X=X, Weight.matrix=genout)
summary(mout)
#Let's determine if balance has actually been obtained on the variables of interest
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)


### 2 ###
## Breakout instruction notes ##
foo <- read.csv("https://course-resources.minerva.kgi.edu/uploaded_files/mke/00086677-3767/peace.csv")
# extract relevant columns (adding 108 for logdead)
foo <- foo[, c(6:8, 11:16, 99, 50, 114, 49, 63, 136, 108, 109, 126, 48, 160, 142, 10)]
# remove 2 rows with missing data (there are better ways to handle missing data)
foo <- foo[c(-19, -47), ]
# check that all missing data is gone...
which(is.na(foo) == TRUE)
# take a peek at the data set (identify the columns)
head(foo)
## new code ##
# the logistic regression model without interaction terms:
glm1 <- glm(pbs2s3~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
              exp + decade + treaty + untype4,
            data = foo,
            family = "binomial")
# the logistic regression model with interaction term from paper
glm2 <- glm(pbs2s3~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
              exp + decade + treaty + untype4 + wardur*untype4 ,
            data = foo,
            family = "binomial")
# the logistic regression model with new interaction term
glm3 <- glm(pbs2s3~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
              exp + decade + treaty + untype4 + logdead*untype4,
            data = foo,
            family = "binomial")
# logistic regression with logcost instead of logdead
glm4 <- glm(pbs2s3~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
              exp + decade + treaty + untype4 + logcost*untype4,
            data = foo,
            family = "binomial")
## getting the marginal effect
# creating a treated and control group where all the variables are average except
# duration wich goes up to 315 as in the paper.
treated1 <- data.frame("wartype" = rep(mean(foo$wartype),315), 
                      "logcost" = rep(mean(foo$logcost),315),
                      "wardur" = c(1:315), #1 through 315
                      "factnum" = rep(mean(foo$factnum),315),
                      "factnum2" = rep(mean(foo$factnum),315),
                      "trnsfcap" = rep(mean(foo$trnsfcap),315),
                      "develop" = rep(mean(foo$develop),315),
                      "exp" = rep(mean(foo$exp),315),
                      "decade" = rep(mean(foo$decade),315),
                      "treaty" = rep(mean(foo$treaty),315),
                      "untype4" = rep(1,315)) # this is the treatment group

control1 <- data.frame("wartype" = rep(mean(foo$wartype),315), 
                      "logcost" = rep(mean(foo$logcost),315),
                      "wardur" = c(1:315), #1 through 315
                      "factnum" = rep(mean(foo$factnum),315),
                      "factnum2" = rep(mean(foo$factnum),315),
                      "trnsfcap" = rep(mean(foo$trnsfcap),315),
                      "develop" = rep(mean(foo$develop),315),
                      "exp" = rep(mean(foo$exp),315),
                      "decade" = rep(mean(foo$decade),315),
                      "treaty" = rep(mean(foo$treaty),315),
                      "untype4" = rep(0,315)) # this is the control group

alternate.treated1 <- data.frame("wartype" = rep(mean(foo$wartype),315), 
                                 "logcost" = seq(0,16,length.out = 315),#iterating through logcost
                                 "wardur" = rep(mean(foo$wardur),315), #wardur is constant
                                 "factnum" = rep(mean(foo$factnum),315),
                                 "factnum2" = rep(mean(foo$factnum),315),
                                 "trnsfcap" = rep(mean(foo$trnsfcap),315),
                                 "develop" = rep(mean(foo$develop),315),
                                 "exp" = rep(mean(foo$exp),315),
                                 "decade" = rep(mean(foo$decade),315),
                                 "treaty" = rep(mean(foo$treaty),315),
                                 "untype4" = rep(1,315)) # this is the treatment group

alternate.control1 <- data.frame("wartype" = rep(mean(foo$wartype),315), 
                                 "logcost" = seq(0,16,length.out = 315), #iterating through logcost
                                 "wardur" = rep(mean(foo$wardur),315), # wardur is constant
                                 "factnum" = rep(mean(foo$factnum),315),
                                 "factnum2" = rep(mean(foo$factnum),315),
                                 "trnsfcap" = rep(mean(foo$trnsfcap),315),
                                 "develop" = rep(mean(foo$develop),315),
                                 "exp" = rep(mean(foo$exp),315),
                                 "decade" = rep(mean(foo$decade),315),
                                 "treaty" = rep(mean(foo$treaty),315),
                                 "untype4" = rep(0,315)) # this is the control group

#adding interaction terms
treated2 <- cbind(treated1, "wardur*untype4" = treated1$wardur*treated1$untype4)
control2 <- cbind(control1, "wardur*untype4" = treated1$wardur*treated1$untype4)
treated3 <- cbind(treated1, "logdead" = rep(mean(foo$logdead),315), "logdead*untype4" = rep(mean(foo$logdead),315)*treated1$untype4)
control3 <- cbind(control1, "logdead" = rep(mean(foo$logdead),315), "logdead*untype4" = rep(mean(foo$logdead),315)*treated1$untype4)
alternate.treated3 <- cbind(alternate.treated1, "logcost*untype4" = alternate.treated1$logcost*alternate.treated1$untype4)
alternate.control3 <- cbind(alternate.control1, "logcost*untype4" = alternate.treated1$logcost*alternate.treated1$untype4)


# calculating marginal effects
marginaleffect1 <- predict(glm1, newdata = treated1, type = "response") - predict(glm1, newdata = control1, type = "response")
marginaleffect2 <- predict(glm2, newdata = treated2, type = "response") - predict(glm2, newdata = control2, type = "response")
marginaleffect3 <- predict(glm3, newdata = treated3, type = "response") - predict(glm3, newdata = control3, type = "response")

alternate.marginaleffect1 <-  predict(glm1, newdata = alternate.treated1, type = "response") - predict(glm1, newdata = alternate.control1, type = "response")
alternate.marginaleffect3 <-  predict(glm4, newdata = alternate.treated3, type = "response") - predict(glm4, newdata = alternate.control3, type = "response")

# plotting the original model and model with interaction term from the paper
plot(c(1:315), marginaleffect1, type = "l", ylim = c(0,0.8),xlim = c(0,315),lty = 3,
     xlab = "Duration of wars in months", ylab = "Marginal effects of UN peacekeeping operations")
lines(c(1:315), marginaleffect2, ylim = c(0,0.8), xlim = c(0,315))
text(100,0.2,"Model with interaction term")
text(225,0.47,"Dotted: Original model")

# plotting the original model and model with new interaction term
plot(c(1:315), marginaleffect1, type = "l", ylim = c(0,0.9),xlim = c(0,315),lty = 3,
     xlab = "Duration of wars in months", ylab = "Marginal effects of UN peacekeeping operations")
lines(c(1:315), marginaleffect3, ylim = c(0,0.9), xlim = c(0,315))
text(100,0.7,"Model with new interaction term")
text(225,0.47,"Dotted: Original model")

#plotting with logcost as independent variable instead of wardur
plot(seq(1,16,length.out=315), alternate.marginaleffect1, type = "l", ylim = c(0,0.8),xlim = c(0,16),lty = 3,
     xlab = "logcost of wars", ylab = "Marginal effects of UN peacekeeping operations")
lines(seq(1,16,length.out=315), alternate.marginaleffect3, ylim = c(0,0.8), xlim = c(0,16))
text(7,0.7,"Model with new interaction term")
text(12,0.47,"Dotted: Original model")

### 3 ###
# re-download data to have all variables available
foo <- read.csv("https://course-resources.minerva.kgi.edu/uploaded_files/mke/00086677-3767/peace.csv")
# Original
Tr <- rep(0, length(foo$uncint))
Tr[which(foo$uncint != 0 & foo$uncint != 1)] <- 1

# translated to the structure of the actual data
Tr <- rep(0, length(foo$uncint))
Tr[which(foo$uncint != "None")] <- 1
#add the variable to the dataset
foo$treat <- Tr

### 4 ###
## pre-processing ##
foo <- data.frame("treat" = foo$treat, 
                  "wartype" = foo$wartype,
                  "logcost" = foo$logcost,
                  "wardur" = foo$wardur,
                  "factnum" = foo$factnum,
                  "factnum2" = foo$factnum2,
                  "trnsfcap" = foo$trnsfcap,
                  "exp" = foo$exp,
                  "decade" = foo$decade,
                  "treaty" = foo$treaty,
                  "untype4" = foo$untype4,
                  "pbs2l" = foo$pbs2l,
                  "pbs5l" = foo$pbs5l)
# getting rid of NAs
foo <- foo[c(-4,-16,-19,-47,-84,-93,-98),]
# code pbs2l and pbs5l as 1 and 0 instead of "Success" and "Failure"
effect1 <- rep(0, length(foo$treat))
effect1[which(foo$pbs2l != "Success")] <- 1
foo$pbs2l <- effect1
# since pbs5l includes NAs, we need to make a copy of the dataset without the NAs
effect2 <- rep(0, length(foo$treat))
effect2[which(foo$pbs5l != "Success")] <- 1
foo$pbs5l <- effect2


## Logistic regression: ##
# lenient peace building success after 2 years:
fit2l <- glm(pbs2l~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
              exp + decade + treaty + untype4 ,
            data = foo,
            family = "binomial")
effect2l <- mean(predict(fit2l,newdata = foo[which(foo$treat==1),], type = "response")) - mean(predict(fit2l, newdata = foo[which(foo$treat==0),], type = "response"))
mb2l.1 <- MatchBalance(treat ~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
                         exp + decade + treaty + untype4 ,
                       data = foo,
                       nboots = 500)
effect2l


# lenient peace building success after 5 years:
fit5l <- glm(pbs5l~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
               exp + decade + treaty + untype4,
             data = foo,
             family = "binomial")
effect5l <- mean(predict(fit5l,newdata = foo[which(foo$treat==1),], type = "response")) - mean(predict(fit5l, newdata = foo[which(foo$treat==0),], type = "response"))
mb5l.1 <- MatchBalance(treat ~ wartype + logcost + wardur + factnum + factnum2 + trnsfcap + 
                         exp + decade + treaty + untype4,
                       data = foo,
                       nboots = 500)
effect5l

## p-score matching: ##
# 2 years lenient peace building success 
set.seed(123)
fit1 <- glm(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap + exp + decade + untype4 +
              I(wartype^2) + I(logcost^2) + I(wardur^2) + I(wardur*untype4) + I(factnum*wardur),
            family = binomial,
            data = foo)
mout1 <- Match(Y = foo$pbs2l, X = fit1$fitted, Tr = foo$treat)
mb2l.2 <- MatchBalance(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap  + exp + decade + untype4 +
                         I(wartype^2) + I(logcost^2) + I(wardur^2) + I(wardur*untype4) + I(factnum*wardur),
                       match.out = mout1,
                       data = foo,
                       nboots = 500)
summary(mout1) # no statistically significant result
# bias adjusted:
mout1ba <- Match(Y = foo$pbs2l, X = fit1$fitted, Tr = foo$treat, BiasAdjust = TRUE)
summary(mout1ba) # no statistically significant result


# 5 years lenient peace building success
fit2 <- glm(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap + exp + decade + untype4 +
              I(wartype^2) + I(logcost^2) + I(wardur^2) + I(wardur*untype4) + I(factnum*wardur),
            family = binomial,
            data = foo)
mout2 <- Match(Y = foo$pbs5l, Tr = foo$treat, X = fit2$fitted)
mb5l.2 <- MatchBalance(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap + exp + decade + untype4 +
                         I(wartype^2) + I(logcost^2) + I(wardur^2) + I(wardur*untype4) + I(factnum*wardur),
                       match.out = mout2,
                       data = foo,
                       nboots = 500)
summary(mout2)
# bias-adjusted:
mout2ba <- Match(Y = foo$pbs5l, Tr = foo$treat, X = fit2$fitted, BiasAdjust = TRUE)
summary(mout2ba)


## genetic matching: ##
library(rgenoud)
set.seed(234)
# 2 years lenient peacebuilding success
X2 <- cbind(foo$wartype, foo$logcost, foo$wardur, foo$factnum, foo$factnum2, foo$trnsfcap, foo$exp, foo$decade)
Y2 <- foo$pbs2l
genout2l <- GenMatch(Tr = foo$treat, X = X2, pop.size = 300, max.generations = 30, wait.generations = 2)
mout2l <- Match(Y = Y2, Tr = foo$treat, X = X2, Weight.matrix = genout2l)
mb2l.4 <- MatchBalance(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap + exp + decade,
                       match.out = mout2l,
                       data = foo)
summary(mout2l)
# bias-adjusted: 
mout2lba <- Match(Y = Y2, Tr = foo$treat, X = X2, Weight.matrix = genout2l, BiasAdjust = TRUE)
summary(mout2lba)

# 5 years lenient peacebuilding success
set.seed(345)
X5 <- cbind(foo$wartype, foo$logcost, foo$wardur, foo$factnum, foo$factnum2, foo$trnsfcap, foo$exp, foo$decade)
Y5 <- foo$pbs5l
genout5l <- GenMatch(Tr = foo$treat, X = X5, pop.size = 300, max.generations = 30, wait.generations = 2)
mout5l <- Match(Y = Y5, Tr = foo$treat, X = X5, Weight.matrix = genout5l)
mb2l.4 <- MatchBalance(treat ~  wartype + logcost + wardur + factnum + factnum2 + trnsfcap + exp + decade,
                       match.out = mout5l,
                       data = foo)
summary(mout5l)
# bias-adjusted:
mout5lba <- Match(Y = Y5, Tr = foo$treat, X = X5, Weight.matrix = genout5l, BiasAdjust = TRUE)
summary(mout5lba)