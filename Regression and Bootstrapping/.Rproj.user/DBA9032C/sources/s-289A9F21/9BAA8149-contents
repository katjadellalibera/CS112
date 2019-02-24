## 1 ##
# data generation
# 999 datapoints:
set.seed(1)
data<- data.frame(variable1=runif(999,min=1,max=50))
data$variable2<-50+data$variable1/10+runif(100,min=-2,max=2)
lm_data=lm(variable2~variable1,data=data)
plot(data)
abline(lm_data)
# outlier:
data_with_outlier<- rbind(data,c( variable1=150,variable2=-150))
lm_data_with_outlier=lm(variable2~variable1, data=data_with_outlier)
plot(data_with_outlier)
abline(lm_data_with_outlier,col="green")
abline(lm_data,col="red")
# regression of original data:
summary(lm_data)
# regression 
summary(lm_data_with_outlier)


## 2 ##
library(Matching)
library(arm)
data(lalonde)
library(dplyr)
library(haven)
# selecting the control group
control_group <- lalonde%>%
  filter(treat==0)
# making linear model
lm_lalonde <- glm(re78~age+educ+re74+re75+I(educ*re74)+I(educ*re75)+
                    I(age*re74)+I(age*re75)+I(age*age)+I(re74*re75),data=control_group)
summary(lm_lalonde)
# simulating interval of expected values
set.seed(2)
edu_median<-median(control_group$educ)
re74_median <- median(control_group$re74)
re75_median <- median(control_group$re75)
sim_results <- sim(lm_lalonde,n.sims=10000)
edu_75percent <- quantile(control_group$educ,probs=0.75)
re74_75percent <- quantile(control_group$re74,probs=0.75)
re75_75percent <- quantile(control_group$re75,probs=0.75)

# expected values with variables at median
lower_confidence_1_vector<-rep(0,55-16)
median_confidence_1_vector<-rep(0,55-16)
upper_confidence_1_vector<-rep(0,55-16)
storage_vector<-rep(0,10000)
for (age in 17:55){
    for (i in 1:10000){
      storage_vector[i]<-
      sim_results@coef[i,1]+
      sim_results@coef[i,2]*age+
      sim_results@coef[i,3]*edu_median+
      sim_results@coef[i,4]*re74_median+
      sim_results@coef[i,5]*re75_median+
      sim_results@coef[i,6]*(edu_median*re74_median)+
      sim_results@coef[i,7]*(edu_median*re75_median)+
      sim_results@coef[i,8]*(age+re74_median)+
      sim_results@coef[i,9]*(age*re75_median)+
      sim_results@coef[i,10]*(age*age)+
      sim_results@coef[i,11]*(re74_median*re75_median)
    }
  lower_confidence_1_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[1]
  median_confidence_1_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[2]
  upper_confidence_1_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[3]  
}
lower_confidence_1_vector
upper_confidence_1_vector
confidence_1=data.frame(lower=lower_confidence_1_vector,median=median_confidence_1_vector,
                        upper=upper_confidence_1_vector)
write.csv(confidence_1,file="confidence_1.csv")

# expexted values with variables at 75% quantile
lower_confidence_2_vector<-rep(0,55-16)
median_confidence_2_vector<-rep(0,55-16)
upper_confidence_2_vector<-rep(0,55-16)
storage_vector<-rep(0,10000)
for (age in 17:55){
  for (i in 1:10000){
    storage_vector[i]<-
      sim_results@coef[i,1]+
      sim_results@coef[i,2]*age+
      sim_results@coef[i,3]*edu_75percent+
      sim_results@coef[i,4]*re74_75percent+
      sim_results@coef[i,5]*re75_75percent+
      sim_results@coef[i,6]*(edu_75percent*re74_75percent)+
      sim_results@coef[i,7]*(edu_75percent*re75_75percent)+
      sim_results@coef[i,8]*(age+re74_75percent)+
      sim_results@coef[i,9]*(age*re75_75percent)+
      sim_results@coef[i,10]*(age*age)+
      sim_results@coef[i,11]*(re74_75percent*re75_75percent)
  }
  lower_confidence_2_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[1]
  median_confidence_2_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[2] 
  upper_confidence_2_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5,0.975))[3]  
}
confidence_2=data.frame(lower=lower_confidence_2_vector,median=median_confidence_2_vector
                        ,upper=upper_confidence_2_vector)
write.csv(confidence_2,file="confidence_2.csv")

# plot for expected values
plot(x=c(1:13000),y=c(1:13000),type="n", xlim = c(16,56),ylim = c(1,13000),
     main = "Plot of the expected value of re78 dependent on age",
     xlab= "Age", ylab="Expected value of re78 (95% confidence)")

for (age in 17:55){
  segments(
    x0= age,
    y0=confidence_1$lower[age-16],
    x1=age,
    y1=confidence_1$upper[age-16],
    lwd=2
    )
}

for (age in 17:55){
  segments(
    x0= age,
    y0=confidence_2$lower[age-16],
    x1=age,
    y1=confidence_2$upper[age-16],
    lwd=2
  )
}

# predictions with values constant and sigmas
lower_confidence_3_vector<-rep(0,55-16)
median_confidence_3_vector<-rep(0,55-16)
upper_confidence_3_vector<-rep(0,55-16)
storage_vector<-rep(0,10000)
for (age in 17:55){
  for (i in 1:10000){
    storage_vector[i]<-
      sim_results@coef[i,1]+
      sim_results@coef[i,2]*age+
      sim_results@coef[i,3]*edu_median+
      sim_results@coef[i,4]*re74_median+
      sim_results@coef[i,5]*re75_median+
      sim_results@coef[i,6]*(edu_median*re74_median)+
      sim_results@coef[i,7]*(edu_median*re75_median)+
      sim_results@coef[i,8]*(age+re74_median)+
      sim_results@coef[i,9]*(age*re75_median)+
      sim_results@coef[i,10]*(age*age)+
      sim_results@coef[i,11]*(re74_median*re75_median)+
      rnorm(1,0,sim_results@sigma[i])
  }
  lower_confidence_3_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[1]
  median_confidence_3_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[2]  
  upper_confidence_3_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[3]  
}
confidence_3=data.frame(lower=lower_confidence_3_vector,median=median_confidence_3_vector,
                        upper=upper_confidence_3_vector)
write.csv(confidence_3,file="confidence_3.csv")

# predictions with values at 75% quantile and sigmas
lower_confidence_4_vector<-rep(0,55-16)
median_confidence_4_vector<-rep(0,55-16)
upper_confidence_4_vector<-rep(0,55-16)
storage_vector<-rep(0,10000)
for (age in 17:55){
  for (i in 1:10000){
    storage_vector[i]<-
      sim_results@coef[i,1]+
      sim_results@coef[i,2]*age+
      sim_results@coef[i,3]*edu_75percent+
      sim_results@coef[i,4]*re74_75percent+
      sim_results@coef[i,5]*re75_75percent+
      sim_results@coef[i,6]*(edu_75percent*re74_75percent)+
      sim_results@coef[i,7]*(edu_75percent*re75_75percent)+
      sim_results@coef[i,8]*(age+re74_75percent)+
      sim_results@coef[i,9]*(age*re75_75percent)+
      sim_results@coef[i,10]*(age*age)+
      sim_results@coef[i,11]*(re74_75percent*re75_75percent)+
      rnorm(1,0,sim_results@sigma[i])
  }
  lower_confidence_4_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[1]
  median_confidence_4_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[2]
  upper_confidence_4_vector[age-16]<-quantile(storage_vector,probs=c(0.025,0.5, 0.975))[3]  
}
confidence_4=data.frame(lower=lower_confidence_4_vector,median=median_confidence_4_vector,
                        upper=upper_confidence_4_vector)
write.csv(confidence_4,file="confidence_4.csv")

# plot for prediction
plot(x=c(1:13000),y=c(1:13000),type="n", xlim = c(16,56),ylim = c(-8000,19000),
     main = "Plot of the predicted values of re78 dependent on age",
     xlab= "Age", ylab="predicted values of re78 (95% confidence)")

for (age in 17:55){
  segments(
    x0= age,
    y0=confidence_3$lower[age-16],
    x1=age,
    y1=confidence_3$upper[age-16],
    lwd=2
  )
}

for (age in 17:55){
  segments(
    x0= age,
    y0=confidence_4$lower[age-16],
    x1=age,
    y1=confidence_4$upper[age-16],
    lwd=2
  )
}


## 3 ##
library(datasets)
library(boot)
# disregard treatment 2 and add indicator
plant_growth_filtered<- PlantGrowth%>%
  filter(group!="trt2")%>%
  mutate(indicator=as.numeric(group=="trt1"))
# bootstrap the 95% confidence interval
# auxiliary funtion
boot_fn <- function(data,index) return(coef(lm(weight~indicator,
                                               data=plant_growth_filtered),subset=index))
boot_plant_growth<- boot(plant_growth_filtered,boot_fn,10000)
conf_intervals <- apply(boot_plant_growth$t,2,quantile,c(0.025,0.975))
  
## 4 ##
# R^2 function using the fact that it's correlation squared
r_function<- function(Y,pred_y) cor(Y,pred_y)^2

## 5 ##
# import data
library(foreign)
nsw <- read.dta("nsw.dta")
# predictor model
nsw_control<-nsw%>%
  filter(treat==0)
nsw_treat<- nsw%>%
  filter(treat==1)
treatment_observations<-data.frame(age=nsw_treat$age,education=nsw_treat$education,black=nsw_treat$black,
                         hispanic=nsw_treat$hispanic,married=nsw_treat$married,nodegree=nsw_treat$nodegree,
                         re75=nsw_treat$re75)
control_observations<-data.frame(age=nsw_control$age,education=nsw_control$education,black=nsw_control$black,
                                hispanic=nsw_control$hispanic,married=nsw_control$married,nodegree=nsw_control$nodegree,
                                re75=nsw_control$re75)

lm_nsw <- glm(treat~age+education+black+hispanic+married+nodegree+re75,data=nsw)

treatment_probability <- function(coefs,person){
  logit <- coefs[1] + person[1]*coefs[2] +
    person[2]*coefs[3] +
    person[3]*coefs[4] + 
    person[4]*coefs[5] +
    person[5]*coefs[6] +
    person[6]*coefs[7] +
    person[7]*coefs[8]
  return(logit)
}
sim_glm <- sim(lm_nsw, 1000)
storage_vector <- rep(0,1000)
observation_vector <- rep(0,length(treatment_observations))
control_observation_vector<- rep(0,length(control_observations))

# getting probabilities for treatment group
for (observation in 1:length(treatment_observations)){
  for (i in 1:1000){
    storage_vector[i]<-treatment_probability(sim_glm@coef[i,],treatment_observations[observation,])
  }
  data.frame(observation_vector[observation]<-mean(as.numeric(storage_vector)))
}
# getting probabilities for control group
for (observation in 1:length(control_observations)){
  for (i in 1:1000){
    storage_vector[i]<-treatment_probability(sim_glm@coef[i,],control_observations[observation,])
  }
  data.frame(control_observation_vector[observation]<-mean(as.numeric(storage_vector)))
}

#plotting the histogram
library(ggplot2)
ggplot(data=observation_vector, aes(ovservation_vector))+
  geom_histogram()