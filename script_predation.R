#######################################################################################
################  EFFECT OF PREDATION RISK AND TEMPERATURE ON #########################
############################# LOBSTER METABOLIC RATES 

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(lme4)
library(pracma)
library(lattice)
library(lubridate)
library(nlme)
library(mgcv)
library(geepack)
library (gee)
library (MASS)
library(lsmeans)
library(ggthemes)
library(multcomp)
library(effects)

RMR    <-  read.table(file = "RMR.txt", header = TRUE,dec = ".")
scope  <-  read.table(file = "AEROBIC.txt", header = TRUE,dec = ".")

RMR$ftreat   <- factor(RMR$treat)
RMR$ftemp    <- factor(RMR$temp)
RMR$fID      <- factor(RMR$ID)
RMR$fperiod  <- factor(RMR$period)

scope$ftreat <- factor(scope$treat)
scope$ftemp  <- factor(scope$temp)
scope$fID    <- factor(scope$ID)

RMR$nperiod  <- as.numeric(RMR$period)
scope$nrisk  <- as.numeric(scope$treat)

## Explo RMR -------------------------------------------------------------------
source('utils.R')

respo  <- cbind(RMR=RMR$RMR, scope[,c(5:9)])

cor(respo)
cor.prob(respo)
flattenSquareMatrix(cor.prob(respo))
chart.Correlation(respo)

## Model RMR -------------------------------------------------------------------

library(car)  ## para error tipo II de efectos no-balanceados

RMR.md0 <- lme(RMR~ ftreat:fperiod -1,random=~temp|fID,data=RMR) ## As suggested the other day, we specify errors from ID by temperature
summary(RMR.md0) 
anova(RMR.md0)
Anova(RMR.md0, type=c("II"))

f1 <-  ggplot(data=subset(RMR, temp==20), aes(nperiod,RMR)) + geom_point(size=3) + 
  facet_grid(. ~ ftreat) + stat_smooth( aes(nperiod,RMR), method = "lm") + 
  scale_x_reverse(breaks = round((c(1.0,2.0)),digits = 3),
                  labels = c('Day','Night')) 
print(f1)

RMR.md20 <- lme(RMR~ ftreat * nperiod,random=~1|fID,data=subset(RMR, temp==20))
summary(RMR.md20)
Anova(RMR.md20, type=c("II"))

f2 <- ggplot(data=subset(RMR, temp==23), aes(nperiod,RMR)) + geom_point(size=3) + 
  facet_grid(. ~ ftreat) + stat_smooth( aes(nperiod,RMR), method = "lm") + 
  scale_x_reverse(breaks = round((c(1.0,2.0)),digits = 3),
                  labels = c('Day','Night')) 
print(f2)

RMR.md23 <- lme(RMR~ ftreat * nperiod,random=~1|fID,data=subset(RMR, temp==23))
summary(RMR.md23)
Anova(RMR.md23, type=c("II"))

## Model RMR -------------------------------------------------------------------

MMR.md1<-lme(MMR~ftemp + nrisk:ftemp -1,random=~1|fID,data=scope) 
Anova(MMR.md1)
anova(MMR.md1)

f3 <- ggplot(data=scope, aes(nrisk,MMR)) + geom_point(size=3) + 
  facet_grid(. ~ ftemp) + stat_smooth( aes(nrisk,MMR), method = "lm") + 
  scale_x_continuous(breaks = round((c(1.0,2.0)),digits = 3),
                     labels = c('No Risk','Risk')) 
print(f3)


## Model AMR -------------------------------------------------------------------

AMR.md1<-lme(AMR~nrisk*ftemp,random=~1|fID,data=scope) 
Anova(AMR.md1)
anova(AMR.md1)

f4 <- ggplot(data=scope, aes(nrisk,AMR)) + geom_point(size=3) + 
  facet_grid(. ~ ftemp) + stat_smooth( aes(nrisk,AMR), method = "lm") + 
  scale_x_continuous(breaks = round((c(1.0,2.0)),digits = 3),
                  labels = c('No Risk','Risk')) 
print(f4)


## Model AS -------------------------------------------------------------------

AS.md1<-lme(AS~nrisk*ftemp,random=~1|fID,data=scope) 
Anova(AS.md1)
anova(AS.md1)

f5 <- ggplot(data=scope, aes(nrisk,AS)) + geom_point(size=3) + 
  facet_grid(. ~ ftemp) + stat_smooth( aes(nrisk,AS), method = "lm") + 
  scale_x_continuous(breaks = round((c(1.0,2.0)),digits = 3),
                     labels = c('No Risk','Risk')) 
print(f5)

## Model SP -------------------------------------------------------------------

SP.md1<-lme(SP~ftemp + nrisk:ftemp -1,random=~1|fID,data=scope) 
Anova(SP.md1)
anova(SP.md1)

f6 <- ggplot(data=scope, aes(nrisk,SP)) + geom_point(size=3) + 
  facet_grid(. ~ ftemp) + stat_smooth( aes(nrisk,SP), method = "lm") + 
  scale_x_continuous(breaks = round((c(1.0,2.0)),digits = 3),
                     labels = c('No Risk','Risk')) 
print(f6)
 

























