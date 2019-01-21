library(ggplot2)
library(dplyr)

### Multilateral Development Institution Data
mdid <- read.csv("https://tinyurl.com/yb4phxx8") # read in the data

#columns representing dates
date.columns <- c(11, 12, 14, 15, 16, 17, 18, 25)

#replacing missing values with NA in date columns
for(i in date.columns)
  {
  
  which_values_are_missing <- which(as.character(mdid[, i]) == "")
  mdid[which_values_are_missing, i] <- NA
  mdid[, i] <- as.Date(as.character(mdid[, i]))
}

#determine which rows have no circulation date given
which.have.NAs <- which(is.na(mdid$CirculationDate == TRUE))
#filter missing values and values with dates before 2008-01-01
new_mdid <- mdid[-which.have.NAs, ]%>%
  filter(CirculationDate >= "2008-01-01")

## 1 ##
# a #
#histogram of the time differences
new_mdid$original.duration= new_mdid$OriginalCompletionDate - new_mdid$ApprovalDate
nas<- which(is.na(new_mdid$original.duration == TRUE))
hist(as.numeric(new_mdid$original.duration[-nas]),
     main="Histogram of the time difference between
     original completion date and approval date",
     xlab="time difference in days")

#plot the time difference over the completion date
ggplot(new_mdid, aes(x=CirculationDate, y=original.duration))+
    geom_point()

#other analysis
mean(new_mdid$original.duration[-nas])
median(new_mdid$original.duration[-nas])
quantile(new_mdid$original.duration[-nas])

# b #
#actual duration of the projects
new_mdid$actual.duration=new_mdid$RevisedCompletionDate - new_mdid$ApprovalDate

#plot of actual time differences
ggplot(new_mdid, aes(x=CirculationDate, y=actual.duration))+
  geom_point()

## 2 ##
# getting the rating percentages
round(100*prop.table(table(new_mdid$Rating)),digits=0)

## 3 ##
# exclude PPTA projects
non_PPTA<- new_mdid%>%
  filter(Type!="PPTA")
# calculate rating percentages
round(100*prop.table(table(non_PPTA$Rating)),digits=0)

## 4 ## 
