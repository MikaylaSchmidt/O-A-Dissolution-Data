################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2024 Mikayla Schmidt
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
#Full Thesis Modelling (Ch1+2)
#for looking at allllll the data together

#0.0 ch1 load
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell1 <- read.csv('./Ch1Clean_MasterSet.csv')
dShell1 <- dShell1[,-c(1,15,23)]
dShell1$chapter <- 'ch1'


#ch2 load
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell2 <- read.csv('./Ch2Clean_MasterSet.csv')
dShell2 <- dShell2[,-c(1, 15,16,19,20,23,24,25,28,32)]
dShell2$chapter <- 'ch2'
dShell2$waxYN <- 'No Wax'

#bind into one
masterSet <- rbind(dShell1, dShell2)

#masterset but with ch1expt3 removed
masterFinal <- subset(masterSet, masterSet$exptID != 'T3.1')
masterFinal <- subset(masterFinal, masterFinal$exptID != 'T3.2')
masterFinal <- subset(masterFinal, masterFinal$exptID != 'T3.3')

#fix porosity of less than 1 
#subPor <- which(masterFinal$poros <= 0)
#masterFinal[subPor, 'poros'] <- 0.0001
subPor2 <- which(masterFinal$poros >=1)
masterFinal[subPor2, 'poros'] <- 1


#load in library
library(lmerTest)
library(MuMIn)


#1.0 basic things looking at whether you need to log the data
#yes, need to log data
plot(cMass ~ finalSA, data=masterFinal)
only1 <- lm(cMass ~ finalSA, data=masterFinal)
plot(only1)
summary(only1)
#logged data itself actually pretty okay, but previously proved we need to change to %/% form
plot(log2(cMass) ~ log2(finalSA), data=masterFinal)
only1.1 <- lm(log2(cMass) ~ log2(finalSA), data=masterFinal)
plot(only1.1)
summary(only1.1)

#change this to percentTotal and percentSATotal form 
#add a way to look at individual points, plus a plot of
plot(percentTotal~percentSATotal, data= masterFinal)
plot(log2(percentTotal) ~ log2(percentSATotal), data=masterFinal, xlab ='Logged Standardized Surface Area', ylab = "Logged Standardized Mass Lost")
line1 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=masterFinal)
abline(line1)
summary(line1)
abline(a=0, b=1, lty= 'dashed')


#2.0 cal3 is the kitchen sink one. nothing is significant except for percentSATotal and taxon
#i dont know how to factor in magensium because for the aragonite it's an NA

relMasterFinal <- masterFinal[,c('cMass', 'pMass', 'mass1', "mass2", 'percentTotal', 'finalSA', 'percentSATotal', 'volume', 'cSize', 'poros', 'deviation', 'taxon', 'polymorph', 'exptID', 'pH.01', 'pH.99', 'meanTemp', 'deltaPH', 'duration', 'phHour')]

#quick check for correlation
allNum1 <-masterFinal[,c('percentSATotal', 'volume', 'cSize', 'poros', 'deviation', 'pH.01', 'meanTemp', 'deltaPH', 'duration')]
cor1 <- cor(allNum1)
print(cor1)
xtable(cor1)
cor1 <- cor1[(cor1 > 0.70)]
print(cor1)

cal3 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp)+ log2(pH.01) + log2(duration) + (1|polymorph) + (1|exptID) + (1|taxon), data = relMasterFinal)
  summary(cal3)
  rand(cal3)
  r.squaredGLMM(cal3)
  coef(cal3)
  
  #beta coefficients
  masterScale <- relMasterFinal
  masterScale[,c(5:11,15:20)] <- scale(log2(masterScale[,c(5:11,15:20)]))
  
  cal3beta <- lmer(percentTotal ~  percentSATotal +  volume + poros + deviation + meanTemp + pH.01 + duration + (1|polymorph) + (1|taxon) + (1|exptID), data = masterScale)
  summary(cal3beta)
  round((summary(cal3beta)$coefficients)^2, digits=5)
  
  
#full model minus taxon or expt ID
  cal3.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(pH.01) + log2(duration)  + (1|polymorph), data = masterFinal)
  summary(cal3.1)
  rand(cal3.1)
  r.squaredGLMM(cal3.1)

#best model overall for the full thesis, no taxon
  calBest <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(poros) + log2(deviation) + (1|polymorph), data = masterFinal)
  summary(calBest)
  rand(calBest)
  r.squaredGLMM(calBest)
  coef(calBest)
  
  
  cal3.1beta <- lmer(percentTotal ~  percentSATotal + poros + deviation + (1|polymorph) , data = masterScale)
  summary(cal3.1beta)
  round((summary(cal3.1beta)$coefficients)^2, digits=3)

  
  cal3.1Best <- lmer(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation)  + (1|polymorph), data = masterFinal)
  summary(calBest)
  rand(calBest)
  r.squaredGLMM(calBest)
  
 #how represent this visually: plot by polymorph
boxplot(perMassSA ~ taxon, data=masterFinal)  
boxplot(perMassSA ~ taxon, data=masterFinal, cex.axis=0.75)
points(massSA ~ taxon, data=masterFinal)  

#plot(lm((log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation) + polymorph, data = masterFinal)))

plot(log2(percentTotal) ~  log2(percentSATotal) * polymorph, data = masterFinal)

#boxplot of the full dataset - not really very nice, shows a lot of error, very smushed
boxplot(perMassSA ~ polymorph, data=masterFinal, cex.axis=0.75, xlab = "Polymorph", ylab= "Standardized Mass Lost/Standardized Surface Area", ylim=c(0,2.5)) 
stripchart(perMassSA ~ polymorph, data = masterFinal, vertical=TRUE, add=TRUE, pch=1)


#plot time
pdf('fullSA.pdf', width =7, height=7)
plot(log2(percentTotal) ~ log2(percentSATotal), data = masterFinal, type='n', ylab = 'Log2 (Standardized Mass Lost)', xlab='Log2 (Standardized Surface Area)',ylim=c(-1.2,3.2),xlim=c(-1.2,3.2) )

masterFinalArag <- subset(masterFinal, masterFinal$polymorph == 'Aragonite')
points(log2(percentTotal) ~ log2(percentSATotal), data = masterFinalArag, pch=1, col='blue')
line1 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=masterFinalArag)
abline(line1, col='blue', lwd=2)
summary(line1)

masterFinalCalc <- subset(masterFinal, masterFinal$polymorph == 'Calcite')
points(log2(percentTotal) ~ log2(percentSATotal), data = masterFinalCalc, pch=2, col='red')
line2 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=masterFinalCalc)
abline(line2, col='red', lwd=2)
summary(line2)

masterFinalMixed <- subset(masterFinal, masterFinal$polymorph == 'Mixed')
points(log2(percentTotal) ~ log2(percentSATotal), data = masterFinalMixed, pch=0, col= '#499E76')
line3 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=masterFinalMixed)
abline(line3, col='#499E76', lwd=2)
summary(line3)

abline(a=0,b=1, col='black', lwd=2)

legend(x= 'topleft', legend = c('Aragonite', 'Calcite', 'Mixed', 'Full Model'), col = c("blue","red",'#499E76','black'), pch=c(1,2,0,NA),lty=1, lwd=2 )  
 dev.off() 
  
    
#3.0  only exptID.
#exptID only has 1% difference on r2 and soaks up the effect of meanTemp
#this also shows how variable for taxon absorbs stuff from volume, poros, and deviation
  
  #below is with meantemp
  cal3.2 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + (1|exptID), data = masterFinal)
  summary(cal3.2)
  rand(cal3.2)
  r.squaredGLMM(cal3.2)
  
  #and then without expt variables
  cal3.3 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + (1|exptID), data = masterFinal)
  summary(cal3.3)
  rand(cal3.3)
  r.squaredGLMM(cal3.3)

  
#4.0 only taxon, kitchen sink model minus exptID and polymorph, to look at just taxon. taxon very very strong
  cal4 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + (1|taxon), data = masterFinal)
  summary(cal4)
  rand(cal4)
  r.squaredGLMM(cal4)
  
  #then without measured variables other than percentSATotal, meanTemp very slightly positive
  cal4.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + (1|taxon), data = masterFinal)
  summary(cal4.1)
  rand(cal4.1)
  r.squaredGLMM(cal4.1)

  
#5.0 Best model with and without taxon
  #without taxon
  cal5 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation) + (1|polymorph), data = masterFinal)
  summary(cal5)
  rand(cal5)
  r.squaredGLMM(cal5)  
  
  #with taxon
  cal5.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + poros + (1|taxon), data = masterFinal)
  summary(cal5.1)
  rand(cal5.1)
  r.squaredGLMM(cal5.1)  
  

  
  
  
  
######6.0  Now onto John model looking at % of specimen mass lost
  #first, the full kitchen sink model
  cal7 <- lmer(log2(pMass) ~  log2(finalSA) + poros + log2(deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(duration) + (1|exptID) + (1|taxon) + (1|polymorph), data = masterFinal)
  summary(cal7)
  rand(cal7)
  r.squaredGLMM(cal7)
  
  #same but minus taxon
  cal7.1 <- lmer(log2(pMass) ~  log2(finalSA)  + poros + log2(deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(duration) + (1|exptID) + (1|polymorph), data = masterFinal)
  summary(cal7.1)
  rand(cal7.1)
  r.squaredGLMM(cal7.1)

  #same but only polymorph
  cal7.2 <- lmer(log2(pMass) ~  log2(finalSA) + poros + log2(deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(duration) + (1|polymorph), data = masterFinal)
  summary(cal7.2)
  rand(cal7.2)
  r.squaredGLMM(cal7.2)
  
#6.1 removing taxon/polymorph from the equation to assess exptID effect
  testModel5<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|exptID), data=masterFinal)
  summary(testModel5)
  rand(testModel5)
  r.squaredGLMM(testModel5)
  
  #exptID yes, expt variables and taxon/polymorph no
  testModel5.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + (1|exptID), data=masterFinal)
  summary(testModel5.1)
  rand(testModel5.1)
  r.squaredGLMM(testModel5.1) 

  
#6.2 second, removing exptID/polymorph from the equation to assess taxon effect
  testModel6<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|taxon), data=masterFinal)
  summary(testModel6)
  rand(testModel6)
  r.squaredGLMM(testModel6)
  
  #taxon yes, measured variables and exptID no
  testModel6.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|taxon), data=masterFinal)
  summary(testModel6.1)
  rand(testModel6.1)
  r.squaredGLMM(testModel6.1) 
  
#6.3 Final best model without/with taxon
  #without taxon, with exptID
  best1 <- lmer(log2(pMass) ~  log2(finalSA) + poros + log2(deviation) + log2(volume) +(1|exptID) + (1|polymorph), data = masterFinal)
  summary(best1)
  rand(best1)
  ranef(best1)
  r.squaredGLMM(best1)
  
  #without taxon or exptID
  best1.1 <- lmer(log2(pMass) ~  log2(finalSA) + poros + log2(deviation) + log2(volume) + (1|polymorph), data = masterFinal)
  summary(best1.1)
  rand(best1.1)
  r.squaredGLMM(best1.1)
  
  #with taxon
  best2 <- lmer(log2(pMass) ~  log2(finalSA) + poros + log2(volume) + (1|taxon), data = masterFinal)
  summary(best2)
  rand(best2)
  r.squaredGLMM(best2)
  