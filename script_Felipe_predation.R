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
RMR.md0.1 <- gam(RMR ~ fperiod + ftemp, data=RMR)
par(mfrow=c(2,1))
plot(RMR.md0.1, all.terms=T, residuals=TRUE)

RMR.md0<-gls(RMR~ ftreat*fperiod+ftemp,data=RMR)
RMR.md1<-lme(RMR~ ftreat*fperiod+ftemp,random=~1|fID,data=RMR) ## ID contributes with lots of heterogeneity, so I added as random


anova(RMR.md0,RMR.md1) ## Comparing the inclusion of ID as randome effect, improves the model so I will keep going with mixed models

## Model #2: Juan Carlos's suggestion
RMR.md2<-lme(RMR~ ftreat:fperiod - 1,random=~temp|fID,data=RMR) ## As suggested the other day, we specify errors from ID by temperature
summary(RMR.md2) ### Just concerned that all was significant!
anova(RMR.md2)

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
MMR.md1<-lme(MMR~ftreat*ftemp,random=~1|fID,data=scope) ## ID contributes with lots of heterogeneity, so I added as random
anova(MMR.md0,MMR.md1) ## Comparing the inclusion of ID as randome effect, doesn't improve the model though

## Model #2: Juan Carlos's suggestion
MMR.md3<-lm(MMR~ftreat:ftemp - 1,data=scope)
summary(MMR.md3) ### All was significant, smells a bit odd!
anova(MMR.md3)

### Model - AMR
AMR.md0<-gls(AMR~ftreat*ftemp,data=scope)
AMR.md1<-lme(AMR~ftreat*ftemp,random=~1|fID,data=scope) 
anova(AMR.md0,AMR.md1)
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

