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

#No calcite analysis
#0.0 Setup script 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./Ch1Clean_MasterSet.csv')
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

#Creating new data frame that removes all specimens where polymorph= calcite
pData <- dShell
pData <- subset(pData, pData$polymorph == 'Aragonite')


TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=c("#004949","#009292","#ff6db6","#ffb6db","#490092","#6db6ff","#920000","#924900","#db6d00"), tLine = 1:length(taxa))
TAXA2$taxon <- factor(c('Ethalia','Notocochlis','Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria','Halimeda','Aragonite'))
TAXA2$tPoint <- factor(c('E','N','L', 'T', 'A', 'P', 'F','H','A'))
pData<-merge(pData,TAXA2,by='taxon')


#1.1 Combining Expt 1, Expt 2, and Expt 3 and then subsetting by experiment
#Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
Expt123 <- pData
Expt123$exptNo <- 'Expt1'
Expt123[(Expt123$exptID == 'T3.1'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T3.2'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T3.3'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T2.1'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.2'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.3'),'exptNo'] <-'Expt2'

factor(Expt123$exptNo)
Expt123$exptNo <- as.factor(Expt123$exptNo)
str(Expt123)

#summary and residuals plotting
summary(lm(percentTotal~percentSATotal, data = Expt123))
plot(lm(percentTotal~percentSATotal, data = Expt123))

#needs to be logged
pdf('firstTest.pdf', width =8, height=7)
par(mfrow = c(2,2))
summary(lm(log2(percentTotal)~log2(percentSATotal), data = Expt123))
plot(lm(log2(percentTotal)~log2(percentSATotal), data = Expt123))
dev.off()

#cmass logged?
summary(lm(log2(cMass)~log2(finalSA), data = Expt123))
plot(lm(log2(cMass)~log2(percentSATotal), data = Expt123))
plot(log2(cMass)~log2(percentSATotal), data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2))

#1.2 visual plot of the logged data

pdf('ch1onlySA.pdf', width = 7, height=7)
plot(log2(percentTotal)~log2(percentSATotal), data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = 'Log2 (Percent of Total Surface Area)', ylab = "Log2 (Percent of Total Mass Lost)", ylim=c(-1,3.15), xlim=c(-1,3.15))
line <-lm(log2(percentTotal)~log2(percentSATotal), data = Expt123)
abline(line, lwd=1.75)
#legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
abline(a=0, b=1, lty='dashed', lwd=1.75)
legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= paste(TAXA2$tPoint), text.col= TAXA2$tColor)


dev.off()
#ablines for all of the individual taxa
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- Expt123[(Expt123$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(log2(percentTotal) ~ log2(percentSATotal), data=temp)
  p <- p + 1
  summaryData[p,'r2'] <- round(summary(lm.temp)$adj.r.squared,3)
  summaryData[p, 'p'] <- round(summary(lm.temp)$coefficients[2,4],7)
  summaryData[p,'slope'] <- round(summary(lm.temp)$coefficients[2,1],3)
  summaryData[p,'slope.err'] <- round(summary(lm.temp)$coefficients[2,2],3)
  summaryData[p,'b.est'] <- round(summary(lm.temp)$coefficients[1,1],3)
  summaryData[p,'b.err'] <- round(summary(lm.temp)$coefficients[1,2], 3)
  summaryData[p,'n'] <- nrow(temp)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topleft', legend= paste(TAXA2$taxon,'slope =',summaryData$slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

dev.off()
#2.0 Modelling predictors
#default plot
model0 <- lm(log2(percentTotal) ~ log2(percentSATotal), data = Expt123)
summary(model0)
plot(model0)

Expt123$resid <- residuals(model0)
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of lm(percentTotal ~ percentSATotal)')
abline(h=0)

#create new data frame using only the relevant data aka relData
relData <- Expt123[,c('cMass', 'pMass', 'mass1', "mass2", 'percentTotal', 'finalSA', 'percentSATotal', 'volume', 'cSize', 'thick', 'poros', 'deviation', 'taxon', 'exptID', 'pH.01', 'pH.99', 'meanTemp', 'deltaPH', 'duration', 'phHour')]
noTaxon <- relData[,c('cMass', 'pMass', 'mass1', "mass2", 'percentTotal', 'finalSA', 'percentSATotal', 'volume', 'cSize', 'thick', 'poros', 'deviation', 'exptID', 'pH.01', 'pH.99', 'meanTemp', 'deltaPH', 'duration', 'phHour')]

#correlation check for variables
corDataNum <- Expt123[,c('cMass', 'pMass', 'mass1', "mass2", 'percentTotal', 'finalSA', 'percentSATotal', 'volume', 'cSize', 'thick', 'poros', 'deviation', 'pH.01', 'pH.99', 'meanTemp', 'deltaPH', 'duration', 'phHour')]
cor(corDataNum)

#quick table looking at effect of expt ID
calTry <- lmer(cMass ~  finalSA  + volume + poros + deviation + meanTemp + duration + pH.01 + thick +(1|exptID)+ (1|taxon), data = onlyArag)
summary(calTry)
xtable(round(summary(calTry)$coefficients, digits=3))
rand(calTry)
r.squaredGLMM(calTry)


#load in library
library(lmerTest)
library(MuMIn)
library(cAIC4)

onlyArag <- subset(Expt123, Expt123$polymorph == 'Aragonite')
droplevels(as.factor(onlyArag$polymorph))
droplevels(as.factor(onlyArag$taxon))
onlyArag$taxon <- factor(onlyArag$taxon, levels = c('Alaona',  'Aragonite', 'Ethalia', 'Fustiaria', 'Halimeda', 'Liloa', 'Notocochlis', 'Pinguitellina', 'Turbo'))

#is mass/SA plot better fit by only calcite? 
plot(cMass ~ finalSA, data=onlyArag)
only1 <- lm(cMass ~ finalSA, data=onlyArag)
summary(only1)
plot(log2(cMass) ~ log2(finalSA), data=onlyArag)
only1.1 <- lm(log2(cMass) ~ log2(finalSA), data=onlyArag)
summary(only1.1)

#need to change this to percentTotal and percentSATotal form 
plot(percentTotal~percentSATotal, data= onlyArag)
plot(log2(percentTotal) ~ log2(percentSATotal), data=onlyArag, col=tColor, pch=tPoint)
line1 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=onlyArag)
abline(line1)
summary(line1)

#perMassSA plot fit with just aragonite
#this is bad - if you look at plot of linear model, data is skewed
#so you should log(2) it. see only3 plot, is much better than only2
plot(percentTotal ~ percentSATotal, data=onlyArag, ylab = 'Standardized Mass', xlab='Standardized Surface Area')
only2 <- lm(percentTotal ~ percentSATotal, data=onlyArag)
summary(only2)
plot(only2)

only3 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=onlyArag)
summary(only3)
plot(only3)


#3.0 First, the model looking at percentTotal ~ percentSATotal

  cal1 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(thick) + (1|exptID) + (1|taxon) , data = onlyArag)
  summary(cal1)
  r.squaredGLMM(cal1)
  rand(cal1)
  
  cal2 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(thick) + (1|exptID), data = onlyArag)
  summary(cal2)
  r.squaredGLMM(cal2)
  
  cal2.1 <- lmer(log2(cMass) ~  log2(finalSA) + log2(volume) + poros + log2(deviation) + log2(thick) + (1|exptID), data = onlyArag)
  summary(cal2.1)
  r.squaredGLMM(cal2.1)
  
  #quick correlation matrix for the data
  allNum1 <-onlyArag[,c('percentSATotal', 'volume', 'cSize', 'poros', 'deviation', 'thick', 'pH.01', 'meanTemp', 'deltaPH', 'duration')]
  cor1 <- cor(allNum1)
  print(cor1)
  xtable(cor1)
  cor1 <- cor1[(cor1 > 0.80)]
  print(cor1)
  
  
  #cal3 is the kitchen sink one. nothing is significant except for percentSATotal and taxon. knocked volume and duration due to collinearity
  #cal3 is one of the important figures in ch1
  cal3 <- lmer(log2(percentTotal) ~  log2(percentSATotal)  + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(duration) + log2(pH.01) + log2(thick) +(1|exptID)+ (1|taxon), data = onlyArag)
  summary(cal3)
  rand(cal3)
  r.squaredGLMM(cal3)

  #beta coefficient arag
  aragScale <- onlyArag
  aragScale[,c(15,17:24,27,30:38)] <- scale(log2(aragScale[,c(15,17:24,27,30:38)]))
  
  arag3beta <- lmer(percentTotal ~  percentSATotal +  volume + poros + deviation + meanTemp + pH.01 + duration +thick + (1|exptID) + (1|taxon), data = aragScale)
  summary(arag3beta)
  round((summary(arag3beta)$coefficients)^2, digits=5)
  
  cal3.0 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(thick) + (1|exptID) + (1|taxon), data = onlyArag)
  summary(cal3.0)
  rand(cal3.0)
  r.squaredGLMM(cal3.0)

  #this is your best model for everything
  cal3Best <- lm(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(thick), data = onlyArag)
  dredge(cal3Best)
  summary(cal3.0)
  cal3.0Best <- lm(log2(percentTotal) ~  log2(percentSATotal) + log2(poros) + log2(thick), data = onlyArag)
  summary(cal3.0Best)
  plot(cal3.0Best)
  
  arag3.0Bestbeta <- lm(percentTotal ~  percentSATotal + poros + thick, data = aragScale)
  summary(arag3.0Bestbeta)
  
#3.1 exptID only has 1% difference on r2 and soaks up the effect of meanTemp
#this also shows how variable for taxon absorbs stuff from volume, poros, and deviation
  
  #below is with meantemp
  cal3.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(thick) + (1|exptID), data = onlyArag)
  summary(cal3.1)
  rand(cal3.1)
  r.squaredGLMM(cal3.1)
  
  #and then without expt variables
  cal3.2 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(thick) + (1|exptID), data = onlyArag)
  summary(cal3.2)
  rand(cal3.2)
  r.squaredGLMM(cal3.2)



#3.2 kitchen sink model minus exptID, to look at just taxon. taxon very very strong
  cal4 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(thick) + (1|taxon), data = onlyArag)
  summary(cal4)
  rand(cal4)
  r.squaredGLMM(cal4)
  
  #then without measured variables other than percentSATotal, meanTemp very slightly positive
  cal4.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + (1|taxon), data = onlyArag)
  summary(cal4.1)
  rand(cal4.1)
  r.squaredGLMM(cal4.1)


#3.3 final pared-down model of things
  #without taxon
  cal5 <- lm(log2(percentTotal) ~  log2(percentSATotal), data = onlyArag)
  summary(cal5)
  
  cal5Best <- lm(log2(percentTotal) ~  log2(percentSATotal)+ log2(poros) + log2(deviation) + log2(thick), data = onlyArag)
  summary(cal5Best)
  dredge(cal5Best)
  
  #with taxon
  cal5.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(poros) + log2(deviation) + log2(thick) + (1|taxon), data = onlyArag)
  summary(cal5.1)
  rand(cal5.1)
  r.squaredGLMM(cal5.1)  

  
  
#4.0 john model looking at % mass lost
  #kitchen sink model, need to add in duration
  cal7 <- lmer(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2 (deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(thick) + log2(duration) + (1|exptID) + (1|taxon), data = onlyArag)
  summary(cal7)
  rand(cal7)
  r.squaredGLMM(cal7)
 
  
#4.1 first, removing taxon from the equation
  testModel5<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(thick)+ log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|exptID), data=onlyArag)
  summary(testModel5)
  rand(testModel5)
  r.squaredGLMM(testModel5)
  
  #exptID yes, expt variables and taxon no
  testModel5.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(thick) + (1|exptID), data=onlyArag)
  summary(testModel5.1)
  rand(testModel5.1)
  r.squaredGLMM(testModel5.1)

#4.2 second, removing exptID from the equation 
  testModel6<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(thick)+ log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|taxon), data=onlyArag)
  summary(testModel6)
  rand(testModel6)
  r.squaredGLMM(testModel6)
  
  #taxon yes, measured variables and exptID no
  testModel6.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(deltaPH) + log2(meanTemp) + log2(pH.01) + log2(duration) + (1|taxon), data=onlyArag)
  summary(testModel6.1)
  rand(testModel6.1)
  r.squaredGLMM(testModel6.1)
 
#4.3 
  #best model no taxon
  #removed deltaPH because it only improves fit by 2%
  testModel7 <- lm(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2(deviation) + log2(volume) + log2(thick), data = onlyArag)
  summary(testModel7)
  testModel7.1 <- lm(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2(deviation) + log2(volume) + log2(thick) + log2(deltaPH), data = onlyArag)
  summary(testModel7.1)
  
  #best model with taxon
  #again, can put in deltaPH for an added 3%
  testModel7.2 <- lmer(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2 (deviation) + log2(volume) + log2(thick) + (1|taxon), data = onlyArag)
  summary(testModel7.2)
  rand(testModel7.2)
  r.squaredGLMM(testModel7.2)
  
  testModel7.3 <- lmer(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2(deviation) + log2(volume) + log2(thick) + log2(deltaPH) + (1|taxon), data = onlyArag)
  summary(testModel7.3)
  rand(testModel7.3)
  r.squaredGLMM(testModel7.3)
  