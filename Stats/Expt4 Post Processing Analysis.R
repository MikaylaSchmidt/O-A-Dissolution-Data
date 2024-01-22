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


#Script for analysis of cleaned data after processing
#Can further cut this down by running and looking at plots
#Originally part of script 'Expt1 Analysis2'

setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt4.csv')
source('./PDF Prefs General.R')

#0.0 Measurement checks general
pdf('./measureCheck.pdf', width=6, height=8, page='A4')
par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(3,3,1,1))
TAXA <- unique(pData$taxon)

for (t in TAXA) {
  
  plot(xDim ~ yDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','yDim'))
  plot(xDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','zDim'))
  plot(yDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'yDim','zDim'))
  plot(xDim ~ thick, data=pData[(pData$taxon == t),], main=paste(t,'xDim','thick'))
  plot(xDim ~ mass1, data=pData[(pData$taxon == t),], main=paste(t,'xDim','mass1'))
  
  
}

dev.off()

#1.1 Setup for measurement  by taxa
  pData <- dShell
  TAXA <- unique(pData$taxon)
  pData$taxon <- as.factor(pData$taxon)
  taxa <- sort(unique(pData$taxon))
  tFont <- rep(3,length(taxa))
  tFont[which(taxa == 'Abranda')] <- 1
  taxaAbrev <- substring(taxa,0,5)

#1.2 Dissolution by total mass (mg, cMass) relationship to 4 variables: surface area (mm^2), initial mass, volume, and densityMV
#need to separate taxon out by color rather than letter
pdf('./outFigs/cMassMeasures.pdf', page='A4', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))

#2.0 first, mass lost by taxon 
plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
mtext('Mass lost (mg)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#2.1 Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
plot(cMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
fit <- lm(cMass~finalSA, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)

#2.2 second, plot the mass lost/SA by taxon
if (length(which(!is.na(pData$finalSA))) > 0) {
  plot(cMass/finalSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/finalSA ~ taxon, data=pData)
  mtext('Mass lost / Surface Area (mg/mm\u00b2)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  #   mtext(e)
} else {
  plot(1:1, type='n', ann=FALSE, axes=FALSE)
  mtext('no data', side=1, line=-2)
}

#2.3 plot mass by initial mass
plot(cMass~mass1, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~mass1, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

#2.4 plot volume
plot(cMass~volume, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~volume, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

plot(sqrt(cMass)~sqrt(volume), data=pData)
fit <- glm(sqrt(cMass)~sqrt(volume), data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

#2.5 plot density
plot(cMass~densityMV, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~densityMV, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

#2.6 Mass lost/density by taxon
if (length(which(!is.na(pData$densityMV))) > 0) {
  plot(cMass/densityMV ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/densityMV ~ taxon, data=pData)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  #  mtext(e)
  #following line only needed if you have no data for one part of expt, which you won't
} else {
  plot(1:1, type='n', ann=FALSE, axes=FALSE)
  mtext('no data', side=1, line=-2)
}

dev.off()



#2.7 Modeling this all for cMass - go to Expt4 Stat Analysis after this 

modelcolsC <- c('cMass','taxon', 'finalSA', 'mass1', 'volume', 'densityMV', 'exptID')

fullModelC<-lm(cMass ~.,data=pData[,modelcolsC])
step(fullModelC, direction = 'forward')
step(fullModelC, direction = 'backward')


#2.8 Modeling with CMass, continued
#COMBINE THIS CMASS SECTION WITH THE PREVIOUS SECTION - MANY OF THE SAME PLOTS ARE HERE THAT YOU DO MORE WITH UP THERE, BASICALLY YOU'VE DONE THE SAME THING TWICE.

#1 cMass with finalSA
a5 <- plot(cMass ~ finalSA, data=pData)
a5 <- lm(cMass ~ finalSA, data=pData)
plot(a5)

b5 <- plot(cMass ~ sqrt(finalSA), data=pData)
b5 <- lm(cMass ~ sqrt(finalSA), data=pData)
plot(b5)

pData$rootMass1 <- pData$mass1^(1/3)
c1 <- plot(cMass ~ (rootMass1), data=pData)
c1 <- lm(cMass ~ (rootMass1), data=pData)
plot(c1)

#2 cMass with mass1
a6 <- plot(cMass ~ mass1, data=pData)
a6 <- lm(cMass ~ mass1, data=pData)
plot(a6)

b6 <- plot(log(cMass) ~ log(mass1), data=pData)
b6 <- lm(log(cMass) ~ log(mass1), data=pData)
plot(b6)

#3 cMass with volume
a7 <- plot(cMass ~ volume, data=pData)
a7 <- lm(cMass ~ volume, data=pData)
plot(a7)

b7 <- plot(log(cMass) ~ log(volume), data=pData)
b7 <- lm(log(cMass) ~ log(volume), data=pData)
plot(b7)

#4 cMass with densityMV
a8 <- plot(cMass ~ densityMV, data=pData)
a8 <- lm(cMass ~ densityMV, data=pData)
plot(a8)

b8 <- plot(log(cMass) ~ log(densityMV), data=pData)
b8 <- lm(log(cMass) ~ log(densityMV), data=pData)
plot(b8)



#3.1 pMass: general plots for visualization
#data frame with columns to examine
head(pData)
modelcolsP <- c('pMass','taxon', 'finalSA', 'mass1', 'volume', 'densityMV', 'exptID')

fullModelP<-lm(pMass ~.,data=pData[,modelcolsP])
step(fullModelP, direction = 'forward')
step(fullModelP, direction = 'backward')

#pdf setup
pdf('./outFigs/pMassMeasures.pdf', page='A4', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))

plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='% Mass Lost')
points(pMass ~ taxon, data=pData)
mtext('% Mass Lost', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)

plot(pMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
fit <- glm(pMass~finalSA, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

dev.off() 

#3.2 Trialing statistical analysis
#none of this actually matters, just trying a few things to visual the data
#and see if any of the non linearity issues can be fixed
#1 pMass with finalSA

a1 <- plot(pMass ~ finalSA, data=pData)
a1 <- lm(pMass ~ finalSA, data=pData)
plot(a1)
a1

b1 <- plot(pMass ~ sqrt(finalSA), data=pData)
b1 <- lm(pMass ~ sqrt(finalSA), data=pData)
summary(b1)
plot(b1)

c1 <- plot(sqrt(pMass) ~ sqrt(finalSA), data=pData)
c1 <- lm(sqrt(pMass) ~ sqrt(finalSA), data=pData)
summary(c1)
plot(c1)

pData$finalSAcube <- pData$finalSA^(1/3)
pData$pMasscube <- pData$pMass^(1/3)
d1 <- plot(pMasscube ~ finalSAcube, data=pData)
d1 <- lm(pMasscube ~ finalSAcube, data=pData)
summary(d1)
plot(d1)

e1 <- plot(log(pMass) ~ log(finalSA), data=pData)
e1 <- lm(log(pMass) ~ log(finalSA), data=pData)
summary(e1)
plot(e1)

#2 pMass with mass1

a2 <- plot(pMass ~ mass1, data=pData)
a2 <- lm(pMass ~ mass1, data=pData)
plot(a2)

b2 <- plot(sqrt(pMass) ~ sqrt(mass1), data=pData)
b2 <- lm(sqrt(pMass) ~ sqrt(mass1), data=pData)
summary(b2)
plot(b2)

pData$rootMass1 <- pData$mass1^(1/3)
c2 <- plot(pMass ~ (rootMass1), data=pData)
c2 <- lm(pMass ~ (rootMass1), data=pData)
plot(c2)

d2 <- plot(log(pMass) ~ (rootMass1), data=pData)
d2 <- lm(log(pMass) ~ (rootMass1), data=pData)
plot(d2)

#3 pMass with volume

a3 <- plot(pMass ~ volume, data=pData)
a3 <- lm(pMass ~ volume, data=pData)
plot(a3)

b3 <- plot(sqrt(pMass) ~ sqrt(volume), data=pData)
b3 <- lm(pMass~ sqrt(volume), data=pData)
plot(b3)

#4 pMass with densityMV

a4 <- plot(pMass ~ densityMV, data=pData)
a4 <- lm(pMass ~ densityMV, data=pData)
plot(a4)

b4 <- plot(log(pMass) ~ log(densityMV), data=pData)
b4 <- lm(log(pMass)~ log(densityMV), data=pData)
plot(b4)


