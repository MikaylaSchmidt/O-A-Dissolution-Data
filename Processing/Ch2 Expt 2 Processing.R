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

#Ch 2 Expt 1 Processing
#0.open script, set names and parameters
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 Data")
dShellFull <- read.csv('./Ch2_Expt2.csv', skip=23)


#remove missing or broken taxon
dShell <- subset(dShellFull, exclude == 0)
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

#1.0 surface area calculations based on parallelopiped/rectangular prism area
#this works for all specimens except for urchin spine
#is 2xy + 2yz + 2xz
dShell$finalSA <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#subset of surface area for sealed specs, is just 2xy
subSeal <- which(dShell$seal == 'Y')
dShell[subSeal,'finalSA'] <- 2*(dShell[subSeal,'xDim'] * dShell[subSeal,'yDim'])
  
#And for cone...
#is pi*r(r+sqrt(h^2 +r^2)) 
subCone <- which(dShell$taxon == 'Centrostephanus')
dShell[subCone,'finalSA'] <- pi * dShell[subCone,'yDim']/2 *((dShell[subCone,'yDim']/2)+sqrt((dShell[subCone,'yDim']/2)^2+(dShell[subCone,'xDim'])^2))

#finally, calc for massSA
dShell$massSA <- dShell$cMass/dShell$finalSA

#2.0 volume and density calculations
#Calculations of volume based on assigned shape
#cube volume = lwh 
dShell$volume <- dShell$xDim * dShell$yDim * dShell$zDim

#for cone, volume = pi*r^2*(h/3)
dShell[subCone,'volume'] <- pi * (dShell[subCone,'yDim']/2)^2 * (dShell[subCone,'xDim']/3)

#(mass^1/3)/volume^1/3 with shell volume substituted for shell size 
dShell$densityMV <- dShell$mass1 / dShell$volume

#deviation from spherical and rootMass for an added bonus
dShell$rootMass <- (dShell$mass1)^ (1/3)
dShell$deviation <- (abs(dShell$xDim - dShell$cSize) + abs(dShell$yDim - dShell$cSize) +  abs(dShell$zDim - dShell$cSize))/3
dShell$deviationStan <- dShell$deviation/dShell$cSize


#3.0 calculation of porosity using calcite and 
#pData$poros <- pData$mass1/(pData$volume * 2.83)
#subCal <- which(pData$polymorph =='Calcite')
#pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
#pData$poros <- 1- pData$poros
#need to add final additional bit about anadara: it's mixed calc arag but have to figure out that ratio first


#4.0 export this data as a csv
write.csv(dShell, 'Ch2Clean_Expt2.csv', row.names=TRUE)
dev.off()

