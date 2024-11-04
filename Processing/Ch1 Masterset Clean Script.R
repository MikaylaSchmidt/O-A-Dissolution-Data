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
#Ch1 Masterset Processing Script
#0.Open clean master dataset 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

#remove missing or broken taxon
dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)

#General calculations of mass lost and time
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cTime1 <- strptime(dShell$dateM1, format=dateFormat)
dShell$cTime2 <- strptime(dShell$dateM2, format=dateFormat)
dShell$cDuration <- dShell$cTime2 - dShell$cTime1

#General calculation of shell size for each taxa
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)

#deviation from spherical and rootMass for an added bonus
dShell$rootMass <- (dShell$mass1)^ (1/3)
dShell$deviation <- (abs(dShell$xDim - dShell$cSize) + abs(dShell$yDim - dShell$cSize) +  abs(dShell$zDim - dShell$cSize))/3
dShell$deviationStan <- dShell$deviation/dShell$cSize
pData <- dShell

#natica fix
subNatNew <- which(pData$taxon == 'Notocochlis' & pData$waxYN == 'No Wax')
pData[subNatNew,'finalSA']<- pData[subNatNew, 'finalSA'] * 1.5

#porosity
pData$poros <- pData$mass1/(pData$volume * 2.83)
subCal <- which(pData$polymorph =='Calcite')
pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
pData$poros <- 1- pData$poros

#massSA
pData$massSA <- pData$cMass/pData$finalSA

#%/% form calculation
#%/% calc
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

#time to save
write.csv(pData, 'Ch1Clean_MasterSet.csv', row.names=TRUE)
dev.off()
