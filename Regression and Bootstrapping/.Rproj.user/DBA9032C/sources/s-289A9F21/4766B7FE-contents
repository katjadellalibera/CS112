library(Matching)
data(lalonde)

# validation set approach
library(ISLR)
set.seed(1)
train=sample(392,196)
lm.fit=lm(mpg~horsepower,data=Auto,subset=train)
attach(Auto)
mean((mpg-predict(lm.fit,Auto))[-train]^2)
lm.fit2=lm(mpg~poly(horsepower ,2),data=Auto , subset=train)
mean((mpg -predict (lm.fit2 ,Auto ))[- train]^2)
lm.fit3=lm(mpg~poly(horsepower ,3),data=Auto , subset=train)
mean((mpg -predict (lm.fit3 ,Auto ))[- train]^2)

# LOOCV
glm.fit=glm(mpg~horsepower,data=Auto)
coef(glm.fit)
lm.fit=lm(mpg~horsepower,data = Auto)
coef(lm.fit)

library(boot)
glm.fit=glm(mpg~horsepower,data=Auto)
cv.error=cv.glm(Auto,glm.fit)
cv.error$delta
cv.error=rep(0, 5)
for (i in 1:5){
  glm.fit=glm(mpg~poly(horsepower,i),data=Auto)
  cv.error[i]=cv.glm(Auto,glm.fit)$delta[1]
}
cv.error

#k-fold CV
set.seed(17)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit=glm(mpg~poly(horsepower,i),data=Auto)
  cv.error.10[i]=cv.glm(Auto,glm.fit,K=10)$delta[1]
}
cv.error.10

#bootstrap
alpha.fn=function(data,index){
  X=data$X[index]
  Y=data$Y[index]
  return((var(Y)-cov(X,Y))/(var(X)+var(Y)-2*cov(X,Y)))
}
alpha.fn(Portfolio, 1:100)

set.seed(1)
alpha.fn(Portfolio,sample(100,100,replace=T))

boot(Portfolio,alpha.fn,R=1000)

#lalonde
treatment_group=which(lalonde$treat==1)
control_group=which(lalonde$treat==0)

storage= c()
for (i in 1:1000){
  treatment_group_sample=sample(treatment_group, length(treatment_group),replace=T)
  mean_treatment=mean(lalonde$re78[treatment_group_sample])
  control_group_sample=sample(control_group,length(control_group),replace=T)
  mean_control=mean(lalonde$re78[control_group_sample])
  storage[i]=mean_treatment-mean_control
}
quantile(storage,probs=c(0.025,0.975))
mean(storage)
#poll
quantile(storage, probs=c(0.125,0.875))

