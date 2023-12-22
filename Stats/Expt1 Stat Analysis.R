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

#Statistical Analysis (model selection) using cleaned/processed data
#New Script Data
library(AICcmodavg)
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')

#separate necessary data out from full dataset
#this is for P or C Mass
x <- dShell[c(9,25,26,36,37,38)]

#correlation calculation indicates collinearity
cor(log(x[,c(1,4,5)]))

#Working with only numerical variables' affect on cMass first
#All of the following linear models will be done in log form to make data (roughly) linear
#first, linear model examining final mass (mass1-cMass, aka mass2) by SA, volume, and starting mass
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)))

#then the same but with total mass lost instead of final mass
summary(lm(log(x$cMass)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$cMass)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$cMass)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$cMass)~log(x$mass1)))

#then the same as each of these with additional variable, which density (by mass/volume)
summary(lm(log(x$mass1-x$cMass)~log(x$mass1/x$volume)))
summary(lm(log(x$cMass)~log(x$mass1/x$volume)))

#do it again, only this time including thickness as a variable
x2 <- dShell[c(9,16,25,26,36,37,38)]
summary(lm(log(x2$mass1-x2$cMass)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))
summary(lm(log(x2$cMass)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))

#using mini lso to see which numerical variables are relevant/colinear
minilso<-function(y,X) {
  y <- scale(y)
  X <- scale(X)
  s <- svd(t(X))
  Z <- scale(t(s$u %*% t(s$v)))
  colnames(Z) <- colnames(X)
  l <- lm(y ~ .,data=data.frame(Z))
  return(list(coef = summary(l)$coefficients,vcov = vcov(l)))
}

minilso(log(x2$cMass),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))
minilso(log(x2$mass1 - x2$cMass),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))

#basically, all of this seems to indicate starting mass (mass1), finalSA, thickness, and volume are significant 
#but density (mass/volume) is not significant


#Third data frame with categorical variables for a ....
#Mixed model!
x3 <- data.frame(scale(log(dShell[c(9,16,25,26,36,37,38)])))
summary(lm((x3$cMass)~(x3$mass1)+(x3$finalSA)+(x3$volume)+(x3$thick)))
x3

#linear mixed effects -lmer package
#lmerTest
library(lmerTest)

#first, mixed with taxon
#is significant
summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

#then, mixed with expt ID
summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))

summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
ranova(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
 
#but no R squared until the addition of the MuMIn package
#both taxon and expt ID count for some percentage of the total dissolution
#in total, 79.25% variance explained by the full
library(MuMIn)
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))



#now time for the Percent Mass
#as you can see below, percent mass doesn't work because you can't transform the data to be linear without it losing all meaning
#not normal distribution on histo, even when logged
hist(dShell$pMass)  
hist(log(1-dShell$pMass))

#r squared of mixed model showing differences between models
r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

#summary of full model - correlation is very low
#most likely issue with colinearity
summary(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
