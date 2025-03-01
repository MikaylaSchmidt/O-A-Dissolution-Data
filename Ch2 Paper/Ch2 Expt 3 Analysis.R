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
#0.0 Setup for everything
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell <- read.csv('./Ch2Clean_Expt3.csv')

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Anadara')] <- 1
taxaAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

#1.0 First, cMass relationships  
#Which is mass lost by taxon 
plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
mtext('Mass lost (mg)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
plot(cMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
fit <- lm(cMass~finalSA, data=pData)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)


#2.0 Then, looking at pMass relationships
#percent lost by taxon 
plot(pMass ~ taxon, data=pData[!is.na(pData$pMass),], ylim=c(0,max(pData$pMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(pMass ~ taxon, data=pData[!is.na(pData$pMass),])
mtext('Mass lost (mg)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
plot(pMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (%)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
fit <- lm(pMass~finalSA, data=pData)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)


#3.0 now getting into the actual comparison bit - looking at MassSA
plot(massSA ~ taxon, data=pData)
points(massSA ~ taxon, data=pData)
mtext('Mass lost/surface area (mg/mm2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2) 

#Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
#Same as seen in this script previously but is related to above massSA calc
plot(cMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
fit <- lm(cMass~finalSA, data=pData)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)

#BELOW NEEDS EDITING TO MAKE WORK - give names
#4.0 Power analysis to look at how many more specimens needed of each to tell difference
#matrix of t tests
TAXA3 <- data.frame(taxon=taxa)
ttestResult <- data.frame(taxon = TAXA3$taxon, k = 3, n=NA, f=NA, alpha=0.05, deltaMeanDiff = NA, sigmaSD = NA)
subDefault <- pData

#all of means and things needed for matrix of t tests
expMean <- aggregate(subDefault$massSA, by=list(subDefault$taxon),FUN=mean, na.action=na.omit)
expSD <- aggregate(subDefault$massSA, by=list(subDefault$taxon),FUN=sd)
expSDBig <- max(expSD$x)
expMeanMax <- max(expMean$x)
expMeanMin <- min(expMean$x)
expMeanDiff <- expMeanMax - expMeanMin

#longer script to run power analysis of t test for number of specimens needed to discern massSA differences if present
tax <- vector()
res <- vector()
c <- 0

for (a in 1:nrow(expMean)) {
  for (b in a:nrow(expMean)) {
    
    c <- c + 1
    tax[c] <- paste(expMean[a,'x'],expMean[b,'x'])
    
    if (a == b) {
      res[c] <- -99
    } else {
      
      DELTA <- expMean[a,'x'] - expMean[(b),'x']
      SD <- expSD[a,'x']
      if (expSD[(b),'x'] > SD)
        SD <- expSD[(b),'x']
      
      res[c] <- power.t.test(n=NULL,delta=DELTA,sd=SD, sig.level=0.05,power=0.8)$n
      #res[c] <- power.t.test(n=5,delta=DELTA,sd=SD, sig.level=0.05,power=NULL)$power
    }
  }}

dataTrial <- data.frame(comparison=tax, n=round(res))
dataTrial <- dataTrial[(dataTrial$n > 0),]
dataTrial

##8.0 changing it into %/% standard
cTotal <- aggregate(pData$cMass, by = list(pData$exptID),FUN=sum)
colnames(cTotal)<-c('exptID','cTotal')
pData <- merge(pData, cTotal, by = 'exptID')
pData$percentTotal <- (pData$cMass)/(pData$cTotal) *100

#now for surface area
SATotal <- aggregate(pData$finalSA, by = list(pData$exptID),FUN=sum)
colnames(SATotal)<-c('exptID','SATotal')
pData <- merge(pData, SATotal, by = 'exptID')
pData$percentSATotal <- (pData$finalSA)/(pData$SATotal) * 100

#now standardized %/%
pData$perMassSA <- pData$percentTotal/pData$percentSATotal
plot(perMassSA~taxon, data= pData)
points(perMassSA~taxon, data = pData)  

#plot standard %/%
plot(percentTotal~percentSATotal, data=pData, col= tColor, pch=tPoint)

#quick lm testing
try1 <- lm(cMass~finalSA, data=pData)
summary(try1)
try2 <- lm(percentTotal~percentSATotal, data=pData)
summary(try2)
