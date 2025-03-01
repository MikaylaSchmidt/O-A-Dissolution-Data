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
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)
expt1 <- read.csv('./shellsData_Expt1.csv')
expt1 <- expt1[,-c(1,11,17,18,19,20,21,22,24,25,28,29,30,32,33,34,35,36)]
expt1$waxYN <- 'No Wax'
expt2 <- read.csv('./shellsData_Expt2.csv')
expt2 <- expt2[,-c(1,11,12,20,21,22,26,27,28,29)]
expt3 <- read.csv('./shellsData_Expt3.csv')
expt3 <- expt3[,-c(1,16,17,18,19,20,22,23,26,27,28,30,31,32)]
expt3$waxYN <- 'No Wax'
dShell <- rbind(expt1, expt2, expt3)

#remove missing or broken taxon
dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)

#deviation standardized for cSize
dShell$deviationStan <- dShell$deviation/dShell$cSize
pData <- dShell

#natica fix not needed, already applied to this data
subNatNew <- which(pData$taxon == 'Notocochlis' & pData$waxYN == 'No Wax')
pData[subNatNew,'finalSA']<- pData[subNatNew, 'finalSA'] * 1.5

#porosity
pData$poros <- pData$mass1/(pData$volume * 2.83)
subCal <- which(pData$polymorph =='Calcite')
pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
#pData$poros <- 1- pData$poros

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


#add code for duration, mean temp, 1% pH, 99% pH, and delta pH
#read in experimental data from Hanna logs
hannaData <- read.csv('./Ch 1 ExptData.csv', skip= 20)
hannaData <- hannaData[,-c(2,3,4,6,7,11,12,13,14,15,16,17,18,19)]
pData <- merge(pData, hannaData, by = 'exptID')
pData$deltaPH <- pData$pH.99 - pData$pH.01
#ph over time
pData$phHour <- pData$deltaPH/pData$duration


#add magnesium % for this dataset
pData$mag <- NA
subCal <- which(pData$taxon == 'Calcite')
pData[subCal, 'mag'] <- 0.05
subMargi <- which(pData$taxon == 'Marginopora')
pData[subMargi, 'mag'] <- 0.15


#time to save
write.csv(pData, 'Ch1Clean_MasterSet.csv', row.names=TRUE)
dev.off()

