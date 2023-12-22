################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2023 Mikayla Schmidt
#
#	This program is free software; you can redistribute it and/or modify it 
#	under the terms of version 3 the GNU General Public License as published 
#	by the Free Software Foundation.
#
#	This program is distributed in the hope that it will be useful, but WITHOUT
#	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
#	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
#	more details.
#
#	To view a copy of the license go to:
#	http://www.gnu.org/licenses/licenses.html#GPL
#	
################################################################################


#Script for processing of raw data into full dataset used in analysis
#Originally part of script 'Expt1 Analysis2' and 'Expt 1 Analysis Script 2.0'
#Edited for analysis of Expt 2 data
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt2.csv')
source('./PDF Prefs General.R')
dShell

#setup for individual taxon
#not sure if this is necessary
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
taxaAbrev <- substring(taxa,0,5)

#1.0 look at a few different things to start
plot(cMass~as.factor(taxon), data = dShell)
plot(pMass~as.factor(taxon), data = dShell)

#next one should be all at 0
#something wrong with the -40 and +40 issues.... check those 
#have fixed those, accidentally switched two final masses
plot(mass1~waxMassDiff, data=dShell)

#difference between the cMass with or without wax - should be roughly 1:1
plot(cMass~cMassWax, data=dShell)
abline(a=0, b=1)
a<- lm(cMass~cMassWax, data=dShell)
summary(a)
abline(a)

#now for percentage, which should also be roughly 1:1
plot(pMass~pMassWax, data=dShell)
abline(a=0, b=1)
b<- lm(pMass~pMassWax, data=dShell)
summary(b)
abline(b)

#all of the below only spit out graphs for the shells with wax
plot(as.factor(dShell$taxon), dShell$waxMass1-dShell$waxMass2)
points(as.factor(dShell$taxon), dShell$waxMass1-dShell$waxMass2)

plot(cMassWax~as.factor(taxon), data = dShell)
plot(pMassWax~as.factor(taxon), data = dShell)

#finally, plot of the stiched together version of mass lost
plot(cMassFinal~as.factor(taxon), data = dShell)
plot(pMassFinal~as.factor(taxon), data = dShell)


#2.0 Onto statistical analysis of the actual experiment
library(AICcmodavg)
x <- dShell[c(9,28,29,30,31,32)]

#uh oh! indicated co-linearity
cor(log(x[,c(1,4,5)]))

#Working with only numerical variables' affect on cMass first
#All of the following linear models will be done in log form to make data (roughly) linear
#first, linear model examining final mass (mass1-cMass, aka mass2) by SA, volume, and starting mass
summary(lm(log(x$mass1-x$cMassFinal)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$mass1-x$cMassFinal)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$mass1-x$cMassFinal)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$mass1-x$cMassFinal)~log(x$mass1)))

#then the same but with total mass lost instead of final mass
#from the below, the one following the global model explains the most dissolution, but still only 37%
#this is half as much as in previous expt...
#removing volume doesn't matter, likely because of relationship between SA and volume
summary(lm(log(x$cMassFinal)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$cMassFinal)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$cMassFinal)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$cMassFinal)~log(x$mass1)))

#then the same as each of these with additional variable, which density (by mass/volume)
#density matters for mass final (23.98%) but not total mass lost
summary(lm(log(x$mass1-x$cMassFinal)~log(x$mass1/x$volume)))
summary(lm(log(x$cMassFinal)~log(x$mass1/x$volume)))

#do it again, only this time including thickness as a variable
x2 <- dShell[c(9,17,28,29,30,31,32)]

#again, thickness does matter for final mass (mass2) but not for mass lost 
summary(lm(log(x2$mass1-x2$cMassFinal)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))
summary(lm(log(x2$cMassFinal)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))

#now minilso time
minilso<-function(y,X) {
  y <- scale(y)
  X <- scale(X)
  s <- svd(t(X))
  Z <- scale(t(s$u %*% t(s$v)))
  colnames(Z) <- colnames(X)
  l <- lm(y ~ .,data=data.frame(Z))
  return(list(coef = summary(l)$coefficients,vcov = vcov(l)))
}

minilso(log(x2$cMassFinal),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))
minilso(log(x2$mass1 - x2$cMassFinal),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))

#basically, all of this seems to indicate starting mass (mass1), finalSA, thickness, and volume are significant 
#but density (mass/volume) is not significant

#same thing but scaled
x3 <- data.frame(scale(log(dShell[c(9,17,28,29,30,31,32)])))
summary(lm((x3$cMassFinal)~(x3$mass1)+(x3$finalSA)+(x3$volume)+(x3$thick)))


#3.0 mixed model time
#linear mixed effects -lmer package
#lmerTest
library(lmerTest)

#first, mixed with taxon
#is significant
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

#then, mixed with expt ID
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))

#then with presence of wax
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN)))

#then both expt ID and taxon
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))

#expt ID and waxYN
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN) + (1|dShell$exptID)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN) + (1|dShell$exptID)))

#waxYN and taxon - should redo this after resubsetting taxon
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$waxYN)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$waxYN)))

#all three
summary(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$waxYN) + (1|dShell$exptID)))
rand(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$waxYN) + (1|dShell$exptID)))

#but no R squared until the addition of the MuMIn package
#both taxon and expt ID count for some percentage of the total dissolution
#in total, 63.5% variance explained by the full
library(MuMIn)
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID) + (1|dShell$waxYN)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$waxYN)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$waxYN)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMassFinal~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))
