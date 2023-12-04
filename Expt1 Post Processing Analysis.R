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
#Originally part of script 'Expt1 Analysis2'

setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')
source('./PDF Prefs General.R')

#1.1 Comparison of the Caliper vs ImageJ measurements
#Setup for measurement comp  by taxa
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Aragonite')] <- 1
taxaAbrev <- substring(taxa,0,5)

pdf('./outFigs/ImageJvsCaliper.pdf', page='A4', height = 6, width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0), cex=1)


#1.2 Plot of caliper SA / ImageJ SA with abline at 1
#First plot - anything above can technically all go in prefs script
plot(cSA1/calcSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='', xlim=c(0,11))
points(cSA1/calcSA ~ taxon, data=pData)
abline(a=1, b=0)
mtext('Caliper SA / ImageJ SA', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxaAbrev, font=tFont, cex=0.5, las=2)

#1.3 Log plot of SAs
pData2 <- pData[!is.na(pData$calcSA),]
plot(cSA1~calcSA, data=pData2, pch=substring(pData2$taxon, 0, 2), log='xy', xlab='log(ImageJSA)', ylab='log(CaliperSA)')
abline(a=0, b=1)
levels(pData$taxon)


#1.4 Log plot of SAs by taxon
plot((cSA1-calcSA)/calcSA~taxon, data=pData2, pch=substring(pData$taxon, 0, 2))
log2(((pData2$cSA1-pData2$calcSA)/pData2$calcSA)+1)


dev.off()

#1.5 Mann-Whitney U Tests of each taxon (for ImageJ to caliper measures)
#Subset all 11 taxa for tests

sumTable <- aggregate(dShell$taxon, by=list(dShell$taxon), FUN=length)
sumTable
colnames(sumTable) <- c('taxon', 'n')
sumTable$pValue = 'na'

myWilcox <- function(dShell, taxon, sumTable){
  myRows <- which(dShell$taxon == taxon)
  Abra1<-wilcox.test(dShell[myRows,'cSA1'],dShell[myRows,'calcSA'])
  print(Abra1)
  sumTable[(sumTable$taxon==taxon), 'pValue'] <- round(Abra1$p.value, 7)
  return(sumTable)
}
for(taxon in TAXA)
sumTable <- myWilcox(dShell, taxon=taxon, sumTable)
sumTable

sumTable$pAdjust <- p.adjust(sumTable$pValue, method = 'holm')
write.csv(sumTable, './outTable/sumTable-Wilcox-CaliperVImageJ.csv')

hist(dShell$cMass)

#2.1 Dissolution by total mass (mg, cMass) relationship to surface area (mm^2), initial mass, volume, and densityMV
#need to separate taxon out by color rather than letter

pdf('./outFigs/cMassMeasures.pdf', page='A4', height = 6 , width = 8)
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(4,4,1,1), cex=1)

plot(cMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
fit <- lm(cMass~finalSA, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)

#2.2 second, plot the mass vs calc by taxon

plot(cMass/finalSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
points(cMass/finalSA ~ taxon, data=pData)
mtext('Mass Lost (mg) / Surface Area (mm\u00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxaAbrev, font=tFont, cex=0.5, las=2)

#try <- glm(cMass/finalSA ~ taxon, data=pData)
#summary(try)

#2.3 plot mass by initial mass
plot(cMass~mass1, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~mass1, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)


#2.3 plot volume
plot(cMass~volume, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~volume, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

plot(sqrt(cMass)~sqrt(volume), data=pData)
fit <- glm(sqrt(cMass)~sqrt(volume), data=dShell)
co <- coef(fit)
abline(fit, lwd=2)


#2.4 plot density
plot(cMass~densityMV, data=pData, axes=TRUE, ann=TRUE)
fit <- glm(cMass~densityMV, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

dev.off()


#2.5 Modeling this all for cMass - go to script frogTestShell.R after this 

modelcolsC <- c('cMass','taxon', 'finalSA', 'mass1', 'volume', 'densityMV', 'exptID')

fullModelC<-lm(cMass ~.,data=pData[,modelcolsC])
step(fullModelC, direction = 'forward')
step(fullModelC, direction = 'backward')




#3.1 Second, for pMass
#data frame with columns to examine
head(pData)
modelcolsP <- c('pMass','taxon', 'finalSA', 'mass1', 'volume', 'densityMV', 'exptID')

fullModelP<-lm(pMass ~.,data=pData[,modelcolsP])
step(fullModelP, direction = 'forward')
step(fullModelP, direction = 'backward')

#pdf setup
pdf('./outFigs/pMassMeasures.pdf', page='A4', height = 6 , width = 8)
par(mfrow=c(1,1))

plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='% Mass Lost')
points(pMass ~ taxon, data=pData)
mtext('% Mass Lost', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxaAbrev, font=tFont, cex=0.5, las=2)

plot(pMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
fit <- glm(pMass~finalSA, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)

dev.off() 

#pMass with finalSA
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


#pMass with mass1
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


#pMass with volume
a3 <- plot(pMass ~ volume, data=pData)
a3 <- lm(pMass ~ volume, data=pData)
plot(a3)

b3 <- plot(pMass ~ sqrt(volume), data=pData)
b3 <- lm(pMass~ sqrt(volume), data=pData)
plot(b3)


#pMass with densityMV
a4 <- plot(pMass ~ densityMV, data=pData)
a4 <- lm(pMass ~ densityMV, data=pData)
plot(a4)

b4 <- plot(pMass ~ sqrt(densityMV), data=pData)
b4 <- lm(pMass~ sqrt(densityMV), data=pData)
plot(b4)


#2.5 Then, the same with cMass
#COMBINE THIS CMASS SECTION WITH THE FIRST SECTION - MANY OF THE SAME PLOTS ARE HERE THAT YOU DO MORE WITH UP THERE, BASICALLY YOU'VE DONE THE SAME THING TWICE.


#cMass with finalSA
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

#cMass with mass1
a6 <- plot(cMass ~ mass1, data=pData)
a6 <- lm(cMass ~ mass1, data=pData)
plot(a6)

b6 <- plot(log(cMass) ~ log(mass1), data=pData)
b6 <- lm(log(cMass) ~ log(mass1), data=pData)
plot(b6)


#cMass with volume
a7 <- plot(cMass ~ volume, data=pData)
a7 <- lm(cMass ~ volume, data=pData)
plot(a7)

b7 <- plot(log(cMass) ~ log(volume), data=pData)
b7 <- lm(log(cMass) ~ log(volume), data=pData)
plot(b7)


#cMass with densityMV
a8 <- plot(cMass ~ densityMV, data=pData)
a8 <- lm(cMass ~ densityMV, data=pData)
plot(a8)

b8 <- plot(log(cMass) ~ log(densityMV), data=pData)
b8 <- lm(log(cMass) ~ log(densityMV), data=pData)
plot(b8)




#3.1 Comparison of the pH over time 
#read in data
pH1.1 <- read.csv('./Expt 1.1 Log 1.csv')
pH1.2 <- read.csv('./Expt 1.2 Log 1.csv')
pH1.3 <- read.csv('./Expt 1.3 Log 1.csv')

#convert interval to hours for each (interval # x 5 min / 60)
pH1.1$hours <- (pH1.1$Interval * 5)/60 
pH1.2$hours <- (pH1.2$Interval * 5)/60 
pH1.3$hours <- (pH1.3$Interval * 5)/60 

pdf('./outFigs/pHReps.pdf', page='A4', height = pageHeight, width = pageWidthTwoColumn)
par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(4,4,1,1), cex=1, mgp=c(1.0,0.75,0))

#time ~ pH for replicate 1.1
plot(pH ~ hours, data=pH1.1, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.1 pH Over Time")
fit1.1 <- glm(pH ~ hours, data=pH1.1)
co <- coef(fit1.1)
abline(fit1.1, lwd=2)

#repeat for 1.2
plot(pH ~ hours, data=pH1.2, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.2 pH Over Time")
fit1.2 <- glm(pH ~ hours, data=pH1.2)
co <- coef(fit1.2)
abline(fit1.2, lwd=2)

#and 1.3
plot(pH ~ hours, data=pH1.3, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.3 pH Over Time")
fit1.3 <- glm(pH ~ hours, data=pH1.3)
co <- coef(fit1.3)
abline(fit1.3, lwd=2)

dev.off()

#3.2 Examining pH 
pH1 <- lm(pH ~ hours, data=pH1.1)
summary(pH1)
pH2 <- lm(pH ~ hours, data=pH1.2)
summary(pH2)
pH3 <- lm(pH ~ hours, data=pH1.3)
summary(pH3)
