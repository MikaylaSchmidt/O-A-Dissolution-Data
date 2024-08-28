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

#Ch 2 1.2 Pilot Processing/Power Analysis
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShellFull <- read.csv('./Ch2_Pilot2.csv', skip=15)

#remove missing or broken taxon
dShell <- subset(dShellFull, exclude == 0)

#General calculations of mass lost (total and percent) and cSize
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)

#surface area calculations based on parallelopiped/rectangular prism area
#this works for all pilot specimens except sea urchin spines
#is 2xy + 2yz + 2xz
dShell$finalSA <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#surface area of a cone for urchin spines
#surface area of cone is pi r^2 + pi r (sqrt(h^2 + r^2))
subCone <- which(dShell$shape == 'Cone')
dShell[subCone,'finalSA'] <- (pi * dShell[subCone,'yDim']/2 * dShell[subCone,'zDim']/2) + (pi * dShell[subCone,'yDim']/2 * sqrt((dShell[subCone,'yDim']/2)^2+ dShell[subCone,'xDim']^2))

#now onto the mass
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


#first, subset of data where there is no sealer to look at how the sea urchin spines are going
subDefault <- subset(pData, seal == 'N')

#quick plot of pMass by taxon and massSA by taxon
plot(pMass ~ taxon, data=subDefault)
plot(massSA ~ taxon, data= subDefault)
mean(subDefault$massSA)



#matrix of t tests
TAXA3 <- data.frame(taxon=taxa)
ttestResult <- data.frame(taxon = TAXA3$taxon, k = 3, n=NA, f=NA, alpha=0.05, deltaMeanDiff = NA, sigmaSD = NA)

subCone <- which(ttestResult$taxon == 'Cent'|ttestResult$taxon == 'Ech2')
ttestResult[subCone,'k'] <- 6

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




#2. second, subset focusing on the differences between sealed and unsealed
#minus the sea urchins because they are not sealed
dShell$taxon <- as.factor(dShell$taxon)
taxa <- sort(unique(dShell$taxon))
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
subSeal <- subset(pData, taxon != 'Cent')
subSeal <- subset(subSeal, taxon != 'Ech2')
subSeal$seal <- as.factor(subSeal$seal)
subSeal$taxon <- droplevels(subSeal$taxon)
str(subSeal)

#now time to adjust the surface area based on subseal calcs
subLess <- which((subSeal$taxon == 'Anadara'|subSeal$taxon == 'Oyster'|subSeal$taxon == 'Scallop')  & subSeal$seal == 'Y')
subSeal[subLess,'finalSA'] <- 2 * subSeal[subLess,'xDim'] * subSeal[subLess,'yDim']
subSeal$massSA <- subSeal$cMass/subSeal$finalSA

#boxplot showing the differences  
boxplot(massSA ~ (seal * taxon), data=subSeal)

#now set of t tests to look at whether there is any differences
subAna <- subset(subSeal, subSeal$taxon=='Anadara')
t.test(massSA ~ seal, data=subAna)
subCal <- subset(subSeal, subSeal$taxon=='Calcite')
t.test(massSA ~ seal, data=subCal)
expSD <- aggregate(subSeal$massSA, by=list(subSeal$seal),FUN=sd)
power.t.test(n=3, delta= 0.027 ,sd=0.018, sig.level=0.05, power= NULL, type = c('two.sample'))
subArag <- subset(subSeal, subSeal$taxon=='Aragonite')
t.test(massSA ~ seal, data=subArag)
subOys <- subset(subSeal, subSeal$taxon=='Oyster')
t.test(massSA ~ seal, data=subOys)
subScal <- subset(subSeal, subSeal$taxon=='Scallop')
t.test(massSA ~ seal, data=subScal)

TAXA4 <- data.frame(taxon=taxa)
ttestResult2 <- data.frame(taxon = TAXA4$taxon, k = 5, n=NA, f=NA, alpha=0.05, deltaMeanDiff = NA, sigmaSD = NA)

#now power analysis looking at how many specs you'd need to tell difference
taxa2 <- sort(unique(subSeal$taxon))
TAXA5 <- data.frame(taxon= taxa2)
summaryP2.1 <- data.frame(taxon = TAXA5, n=NA)
p = 0
for (T in TAXA5$taxon) {
  temp <- subSeal[(subSeal$taxon == T),]
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=NULL,delta=abs(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05, power= 0.8, type = c('two.sample'))
  summaryP2.1[p, 'n'] <- power.temp$n
}
summaryP2.1
