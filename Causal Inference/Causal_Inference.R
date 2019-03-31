### 1 ###
## a ##

library(Matching)
data("lalonde")
attach(lalonde)

genout <- GenMatch(Tr=treat, X=X)

summary(mout)

mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75) + re78,
                   match.out=genout, nboots=500)

## b ## 

genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATT", exact = TRUE, M=1,
                   pop.size=16, max.generations=10, wait.generations=1)

#The outcome variable
Y=re78/1000

#
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=treat, X=X, Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#                        
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)

## c ##
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat)

#The outcome variable
Y=re78/1000

#
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=treat, X=X, Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#                        
mb <- MatchBalance(re78~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)


