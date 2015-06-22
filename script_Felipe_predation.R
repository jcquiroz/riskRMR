#######################################################################################
################  EFFECT OF PREDATION RISK AND TEMPERATURE ON #########################
############################# LOBSTER METABOLIC RATES 

#Load packages and library files
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

rm(list=ls(all=TRUE))

#setwd ("C:/Documents and Settings/bricenof/My Documents/analysis_R/Respirometry/predation/JC")
getwd()

#### Analysis #1 - Effect of period (day vs night), predation risk (0-1) and temperature (20-23)
### on routine metabolic rate (RMR). So, period, predation risk and temperature are considered as
### categorical explanatory factors. RMR is oxygen consumed, so continous dependent variable (gamma).

## Note: RMR strongly changes with the circadian pattern resulting in higher RMR during night time compared
##with day time.

## Variables:
### Treatment (treat): ctrl ('no risk') vs pred ('risk')
### ID: lobster individual (eg L1)
### Period: Day (D) vs Night (N)
### Temperature: 20 C vs 23C 

RMR    <-  read.table(file = "RMR.txt", header = TRUE,dec = ".")
scope  <-  read.table(file = "AEROBIC.txt", header = TRUE,dec = ".")


## -----    explo responses  --------------
source('utils.R')

respo  <- cbind(RMR=RMR$RMR, scope[,-c(1:4)])

cor(respo)
cor.prob(respo)
flattenSquareMatrix(cor.prob(respo))
chart.Correlation(respo)

## --------------------------------------------



### Explanatory variables as factors
RMR$ftreat  <- factor(RMR$treat)
RMR$ftemp   <- factor(RMR$temp)
RMR$fID     <- factor(RMR$ID)
RMR$fperiod <- factor(RMR$period)


### Models
### Model #1 
### Note: Interaction is specified, and temperature is only evaluated as level
RMR.md0.1 <- gam(RMR ~ fperiod * ftemp - 1, data=RMR)
anova(RMR.md0.1)
par(mfrow=c(2,1))
plot(RMR.md0.1, all.terms=T, residuals=TRUE)

RMR.md0<-gls(RMR~ ftreat*fperiod+ftemp,data=RMR)
RMR.md01<-glm(RMR~ ftreat*fperiod+ftemp,data=RMR)
RMR.md1<-lme(RMR~ ftreat*fperiod+ftemp,random=~1|fID,data=RMR) ## ID contributes with lots of heterogeneity, so I added as random


anova(RMR.md1) ## Comparing the inclusion of ID as randome effect, improves the model so I will keep going with mixed models

## Model RMR -------------------------------------------------------------------

library(car)  ## para error tipo II de efectos no-balanceados

RMR.md0 <- lme(RMR~ ftreat:fperiod -1,random=~temp|fID,data=RMR) ## As suggested the other day, we specify errors from ID by temperature
summary(RMR.md0) 
anova(RMR.md0)
Anova(RMR.md0, type=c("II"))

f1 <-  ggplot(data=subset(RMR, temp==20), aes(perioo,RMR)) + geom_point(size=3) + 
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


#### Analysis # 2: Effect of predation risk (no risk vs risk) and temperature (20-23) on different metabolic traits:
### Standard metabolic rate (SMR) ## Ojo: This trait wasn't measured under predation risk, so the only effect explored is temperature (SMR~temp)
### Active metabolic rate (AMR)
### Minimun metabolic rate (MMR)
### Aerobic scope rate (AS) : The difference between SMR and AMR
### Aerobic scope rate (AS) : The difference between MMR and AMR

## Explanatory ariables:
### Treatment (treat): ctrl ('no risk') vs pred ('risk')
### Temperature (temp): 20 C vs 23 C
##  ID: lobster individual (eg L1)



### Explanatory variables as factors
scope$ftreat<-factor(scope$treat)
scope$ftemp<-factor(scope$temp)
scope$fID<-factor(scope$ID)

### Model - MMR
MMR.md0<-gls(MMR~ftreat*ftemp,data=scope)
MMR.md1<-lme(MMR~ftreat*ftemp,random=~1|fID,data=scope) 
MMR.md2<-lme(MMR~ ftemp + ftreat:ftemp -1, random=~1|fID,data=scope) ## ID contributes with lots of heterogeneity, so I added as random
anova(MMR.md2) ## Comparing the inclusion of ID as randome effect, doesn't improve the model though
Anova(MMR.md2, type=c("II"))

library(car)
MMR.A <- lm(MMR ~ as.numeric(ftreat):ftemp + ftemp - 1, data=scope) #,
#            contrasts=list(as.numeric(ftreat)=contr.sum, ftemp=contr.sum))
Anova(MMR.A, type=c("II"))


## Model #2: Juan Carlos's suggestion
MMR.md3<-lm(MMR~ftreat:ftemp - 1,data=scope)
summary(MMR.md3) ### All was significant, smells a bit odd!
anova(MMR.md3)

### Model - AMR
AMR.md0<-gls(AMR~ftreat:ftemp -1,data=scope)
AMR.md1<-lme(AMR~ftreat*ftemp,random=~1|fID,data=scope) 
AMR.md1.1<-lme(AMR~ftreat,random=~temp-1|fID,data=scope) 
Anova(AMR.md1)
anova(AMR.md1.1)



### AMR is a single measurement per individual, time isn't involved as other traits. So maybe
## isn't repeated measurements. Lots of heterogenity though

### Model - AS
AS.md0<-gls(AS~ftreat*ftemp,data=scope)
AS.md1<-lme(AS~ftreat*ftemp,random=~1|fID,data=scope) 
anova(AS.md0,AS.md1)
### AS is the difference between AMR (single value) with SMR (minimun values over a specific period, mean)
### ID doesn't improve model!

### Model - SP (Scope for predation)
SP.md0<-gls(SP~ftreat*ftemp,data=scope)
SP.md1<-lme(SP~ftreat*ftemp,random=~1|fID,data=scope) 
anova(SP.md0,SP.md1)
### ID doesn't improve model!

