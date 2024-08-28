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

#Ch 2 1.1 Pilot Processing/Power Analysis
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShellFull <- read.csv('./Ch2_Pilot1.csv', skip=15)

#remove missing or broken taxon
dShell <- subset(dShellFull, exclude == 0)

#General calculations of mass lost (total and percent) and cSize
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)

#surface area calculations based on parallelopiped/rectangular prism area
#this works for all pilot specimens
#is 2xy + 2yz + 2xz
dShell$finalSA <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)
dShell$massSA <- dShell$cMass/dShell$finalSA

#1.0 post cleanup analysis starts here
dShell$taxon <- as.factor(dShell$taxon)
taxa <- sort(unique(dShell$taxon))
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Anadara')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), taxaAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

#quick plot of pMass by taxon and massSA by taxon
plot(pMass ~ taxon, data=pData)
plot(massSA ~ taxon, data= pData)
mean(pData$massSA)

#looking at means
test1 <- lm(pMass~taxon, data=pData)
anova(test1)

test2 <- lm(massSA~taxon, data=pData)
anova(test2)

TAXA3 <- data.frame(taxon=taxa)
ttestResult <- data.frame(taxon = TAXA3$taxon, k = 5, n=NA, f=NA, alpha=0.05, deltaMeanDiff = NA, sigmaSD = NA)


#more about looking at means
expMean <- aggregate(pData$massSA, by=list(pData$taxon),FUN=mean, na.action=na.omit)
expSD <- aggregate(pData$massSA, by=list(pData$taxon),FUN=sd)
expSDBig <- max(expSD$x)
expMeanMax <- max(expMean$x)
expMeanMin <- min(expMean$x)
expMeanDiff <- expMeanMax - expMeanMin
power.anova.test(groups = 5, n = NULL, between.var = expMeanDiff, within.var = expSDBig, sig.level = 0.05, power = 0.80)


#that's just for general anova: telling if there's a difference between any groups. however, between individual groups of specimens...
#have to write a matrix of t tests
expMean <- aggregate(pData$massSA, by=list(dShell$taxon),FUN=mean)
expSD <- aggregate(pData$massSA, by=list(dShell$taxon),FUN=sd)


power.t.test(n=5,delta=(expMean[1,'x']-expMean[2,'x']),sd=expSD[1,'x'], sig.level=0.05,power=NULL)
power.t.test(n=5,delta=(expMean[2,'x']-expMean[3,'x']),sd=expSD[3,'x'], sig.level=0.05,power=NULL)
power.t.test(n=NULL,delta=(expMean[4,'x']-expMean[5,'x']),sd=expSD[5,'x'], sig.level=0.05,power=0.8)

#longer script to run t test for number of specimens needed to discern massSA differences if present
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

#data frame for determining n
dataTrial <- data.frame(comparison=tax, n=round(res))
dataTrial <- dataTrial[(dataTrial$n > 0),]
dataTrial

#data frame for determining power
dataTrial <- data.frame(comparison=tax, power=(res))
dataTrial <- dataTrial[(dataTrial$n > 0),]
dataTrial



#final export as clean data .csv
write.csv(dShell, 'cleanData_Ch2Pilot.csv', row.names=TRUE)
dev.off()