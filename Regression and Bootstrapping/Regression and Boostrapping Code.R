## 1 ##
# data generation
# 999 datapoints:
set.seed(1)
data<- data.frame(variable1=runif(999,min=1,max=50))
data$variable2=50+data$variable1/10+runif(1,min=-20,max=20)
lm_data=lm(variable2~variable1,data=data)
# outlier:
data_with_outlier<- rbind(data,c( variable1=450,variable2=0))
lm_data_with_outlier=lm(variable2~variable1, data=data_with_outlier)
plot(data_with_outlier)
abline(lm_data_with_outlier)
abline(lm_data)
# regression of original data:
summary(lm_data)
# regression 
summary(lm_data_with_outlier)


## 2 ##
library(Matching)
data(lalonde)
library(dplyr)
# selecting the control group
control_group <- lalonde%>%
  filter(treat==0)
# making linear model
lm_lalonde <- glm(re78~age+educ+re74+re75+educ*re74+educ*re75+age*re74+age*re75+age*age+re74*re75,data=control_group)
summary(lm_lalonde)
# simulating interval of expected values


## 3 ##
library(datasets)
# disregard treatment 2 and add indicator
PlantGrowthNew<- PlantGrowth%>%
  filter(group!="trt2")%>%
  mutate(indicator=)
  
## 4 ##
# write R^2 function

## 5 ##
# import data
library(foreign)
nsw <- read.dta("nsw.dta")
# predictor model
