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
### 2 ###
## Breakout instruction notes ##
foo <- read.csv("https://course-resources.minerva.kgi.edu/uploaded_files/mke/00086677-3767/peace.csv")

# extract relevant columns
foo <- foo[, c(6:8, 11:16, 99, 50, 114, 49, 63, 136, 109, 126, 48, 160, 142, 10)]

# remove 2 rows with missing data (there are better ways to handle missing data)
foo <- foo[c(-19, -47), ]

# check that all missing data is gone...
which(is.na(foo) == TRUE)

# take a peek at the data set (identify the columns)
head(foo)


