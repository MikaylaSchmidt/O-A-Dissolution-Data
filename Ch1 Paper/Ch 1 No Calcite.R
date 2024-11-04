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
#Setup script 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

#Creating new data frame that removes all specimens where polymorph= calcite
pData <- dShell
pData <- subset(pData, pData$polymorph == 'Aragonite')


TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3

#Creating new data frame that removes all specimens where polymorph= calcite
pData <- subset(pData, pData$polymorph == 'Aragonite')

#3.1 Combining Expt 1, Expt 2, and Expt 3 and then subsetting by experiment
#Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
Expt123 <- pData
Expt123$exptNo <- 'Expt1'
Expt123[(Expt123$exptID == 'T4.1'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.2'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.3'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T2.1'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.2'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.3'),'exptNo'] <-'Expt2'

factor(Expt123$exptNo)
Expt123$exptNo <- as.factor(Expt123$exptNo)
str(Expt123)

#total carbonate per replicate in mg
cTotal <- aggregate(Expt123$cMass, by = list(Expt123$exptID),FUN=sum)
colnames(cTotal)<-c('exptID','cTotal')
Expt123 <- merge(Expt123, cTotal, by = 'exptID')

#3.2 for plotting
#% calcium carbonate dissolved over total
Expt123$percentTotal <- (Expt123$cMass)/(Expt123$cTotal) *100

#now for surface area
SATotal <- aggregate(Expt123$finalSA, by = list(Expt123$exptID),FUN=sum)
colnames(SATotal)<-c('exptID','SATotal')
Expt123 <- merge(Expt123, SATotal, by = 'exptID')
Expt123$percentSATotal <- (Expt123$finalSA)/(Expt123$SATotal) * 100

#summary and residuals plotting, first OLS then SMA
summary(lm(percentTotal~percentSATotal, data = Expt123))
plot(lm(percentTotal~percentSATotal, data = Expt123))

#visual plot
plot(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentSATotal, data = Expt123)
abline(lm(percentTotal~percentSATotal, data = Expt123))
legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
abline(a=0, b=1, lty='dashed')


#need to update summary data to have new TAXA2 without calcite/margi





#3.5 plot of total predictors
#plot with all of them included EXCEPT volume...
par(mfrow=c(1,1))
model0 <- lm(percentTotal ~ percentSATotal, data = Expt123)
summary(model0)
plot(model0)

Expt123$resid <- residuals(model0)
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of lm(percentTotal ~ percentSATotal)')
abline(h=0)
