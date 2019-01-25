library(ggplot2)
library(dplyr)

### Multilateral Development Institution Data
mdid <- read.csv("https://tinyurl.com/yb4phxx8") # read in the data

#columns representing dates
date.columns <- c(11, 12, 14, 15, 16, 17, 18, 25)

#replacing missing values with NA in date columns
for(i in date.columns)
  {
  
  which.values.are.missing <- which(as.character(mdid[, i]) == "")
  mdid[which.values.are.missing, i] <- NA
  mdid[, i] <- as.Date(as.character(mdid[, i]))
}

#determine which rows have no circulation date given
which.have.NAs <- which(is.na(mdid$CirculationDate == TRUE))
#filter missing values and values with dates before 2008-01-01
new.mdid <- mdid[-which.have.NAs, ]%>%
  filter(CirculationDate >= "2008-01-01")

## 1 ##
# a #
new.mdid$original.duration= difftime(new.mdid$OriginalCompletionDate,
                                    new.mdid$ApprovalDate,
                                    unit = "weeks")
nas<- which(is.na(new.mdid$original.duration == TRUE))
#analyzing common metrics
summary(as.numeric(new.mdid$original.duration[-nas])/52.14)
quantile(as.numeric(new.mdid$original.duration[-nas])/52.14)
#histogram of the time differences
hist(as.numeric(new.mdid$original.duration[-nas])/52.14,
     main= "Histogram of the time difference between
     original completion date and approval date",
     xlab="time difference in years")   


#fit a linear regression on CompletionDate vs duration
fit<- lm(as.numeric(original.duration)/52.14~CirculationDate,data=new.mdid)
summary(fit)
#plot the time difference over the completion date
ggplot(new.mdid, aes(x=CirculationDate, y=as.numeric(original.duration)/52.14))+
  geom_point()+
  xlab("Circulation Date")+
  ylab("Original duation in years")+
  stat_smooth(method="lm",col="red")+
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))


# b #
#actual duration of the projects
new.mdid$actual.duration=difftime(new.mdid$RevisedCompletionDate,
                                  new.mdid$ApprovalDate,
                                  unit="weeks")
nafree.new.mdid<- new.mdid[-nas,]
summary(as.numeric(nafree.new.mdid$actual.duration)/52.14)
#histogram of the acutal durations
hist(as.numeric(nafree.new.mdid$actual.duration)/52.14,
    main= "Histogram of the actual duration",
    xlab="time difference in years") 
#plot both durations on the same plot
ggplot()+
  geom_point(data=nafree.new.mdid, 
             aes(x=CirculationDate, y=as.numeric(original.duration)/52.14,
                color="original duration"))+
  geom_point(data=nafree.new.mdid,
             aes(x=CirculationDate,y=as.numeric(actual.duration)/52.14,
                 color="actual duration"))+
  xlab("Circulation Date")+
  ylab("duration in years")+
  labs(color="type")

#bonus: actual duration vs original duration
fit2<- lm((as.numeric(actual.duration)/52.14) ~ I(as.numeric(original.duration)/52.14),data=new.mdid)
ggplot(new.mdid, aes(x=as.numeric(original.duration)/52.14,y=as.numeric(actual.duration)/52.14))+
  geom_point()+
  xlab("original duration in years")+
  ylab("actual duation in years")

## 2 ##
# getting the rating percentages
percentages<-round(100*prop.table(table(new.mdid$Rating)),digits=0)
percentages
summary(new.mdid$Rating)
pie(percentages,label=c(0,1,2,3))

## 3 ##
# exclude PPTA projects
non.PPTA<- new.mdid%>%
  filter(Type!="PPTA")
# calculate rating percentages
non.PPTA.percentages<-round(100*prop.table(table(non.PPTA$Rating)),digits=0)
non.PPTA.percentages
summary(non.PPTA$Rating)
#make a pie chart
pie(non.PPTA.percentages,label=c(0,1,2,3))

## 4 ## 
#selecting the bottom 25% and top 25% by RevisedAmount
bottom.25<- new.mdid %>%
  arrange(RevisedAmount)%>%
  filter(row_number()/n()<=.25)
top.25<- new.mdid %>%
  arrange(desc(RevisedAmount))%>%
  filter(row_number()/n()<=.25)
#comparing the ratings of the two groups
summary(bottom.25$Rating)
summary(top.25$Rating)
