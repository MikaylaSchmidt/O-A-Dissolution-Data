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
#means file
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
pData <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3
pData$poros <- pData$mass1/(pData$volume * 2.83)
subCal <- which(pData$polymorph =='Calcite')
pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
Expt123 <- pData

abranda <- subset(Expt123, Expt123$taxon =='Alaona')
summary(abranda)
range(abranda$cSize)

#quick test to make sure fig gen script working
#modeltest <-lm(percentTotal~percentSATotal, data = abranda)
#summary(modeltest)


calcite <- subset(Expt123, Expt123$taxon =='Calcite')
summary(calcite)

eth <- subset(Expt123, Expt123$taxon =='Ethalia')
summary(eth)

hali <- subset(Expt123, Expt123$taxon =='Halimeda')
summary(hali)
range(hali$cSize)

liloa <- subset(Expt123, Expt123$taxon =='Liloa')
summary(liloa)

margi <- subset(Expt123, Expt123$taxon =='Marginopora')
summary(margi)
range(margi$cSize)

natica <- subset(Expt123, Expt123$taxon =='Notocochlis')
summary(natica)

pingui <- subset(Expt123, Expt123$taxon =='Pinguitellina')
summary(pingui)

scaph <- subset(Expt123, Expt123$taxon =='Fustiaria')
summary(scaph)

tridacna <- subset(Expt123, Expt123$taxon =='Aragonite')
summary(tridacna)
range(tridacna$cSize)

turbo <- subset(Expt123, Expt123$taxon =='Turbo')
summary(turbo)
range(turbo$cSize)
