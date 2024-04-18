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

#Processing Script from Expt 1 edited for Expt 2
#Script for processing of raw data into full dataset used in analysis
#Originally part of script 'Expt1 Analysis2'


#0.open script, set names and parameters
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShellFull <- read.csv('./shells_Expt2.csv', skip=24)


#remove missing or broken taxon
dShell <- subset(dShellFull, exclude == 0)

dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)

#General calculations of mass lost and time
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)
dShell$cMassWax <- -1 * (dShell$waxMass2 - dShell$waxMass1)
dShell$waxMassDiff <- dShell$cMassWax-dShell$cMass
dShell$pMassWax <- dShell$cMassWax / dShell$mass1
dShell$wax <- dShell$waxMass1 - dShell$mass1


#1.1 Surface area calculations for all of the taxa using caliper measurements
#by default, all are classified as cylinder
dShell$shape <- 'cylinder'
dShell[(dShell$taxon == 'Calcite'),'shape'] <- 'cube'

#pingui and abra technically two domes
dShell[(dShell$taxon == 'Pinguitellina'), 'shape'] <- '2dome'
dShell[(dShell$taxon == 'Abranda'), 'shape'] <- '2dome'

#finally, nat and eth are spherical
dShell[(dShell$taxon == 'Natica'), 'shape'] <-'sphere'
dShell[(dShell$taxon == 'Ethalia'), 'shape'] <-'sphere'



#1.2 surface area caliper calculations
#basic rectangular cube surface area. does not work for halimeda based off how measurements were taken.
#is 2xy + 2yz + 2xz
dShell$finalSA <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#next lines for surface area of domed specimens - abranda and pingui
# is 2 * pi(r^2+h^2)
subDome <- which(dShell$shape == '2dome')
dShell[subDome,'finalSA'] <- 2 * pi * ((dShell[subDome,'xDim']/2 * dShell[subDome,'yDim']/2) + dShell[subDome,'zDim']^2) 

#finally, nat and eth spherical calc
#is 4 * pi * [(x/2 + y/2 + z/2)/3]^2
subSphere <- which(dShell$shape == 'sphere')
dShell[subSphere,'finalSA'] <- 4 * pi * ((dShell[subSphere,'xDim']/2+dShell[subSphere,'yDim']/2 + dShell[subSphere,'zDim']/2)/3)^2


#divide SA by 2 but only for 2domed specimens
subWax <- which(dShell$waxYN == 'Wax' & dShell$shape == '2dome')
dShell[subWax,'finalSA'] <- (dShell[subWax,'finalSA']/2)

#Additionally, combined pMass and cMass stitched together
dShell[subWax, 'cMass'] <- dShell[subWax, 'cMassWax'] 
dShell[subWax, 'pMass'] <- dShell[subWax, 'pMassWax'] 



#1.4 Density as given in Kosnik et al 2009

#Calculations of volume based on assigned shape
#cube volume = lwh 
dShell$volume <- dShell$xDim * dShell$yDim * dShell$zDim
#2dome = dome1 full volume = 1/6 pi h * (3r^2-h^2) 
dShell[subDome,'volume'] <- 1/6 * pi * dShell[subDome,'zDim'] * ((3 * (dShell[subDome,'xDim']/2 * dShell[subDome,'yDim']/2)) + dShell[subDome,'zDim']^2) 
#sphere (nat and eth) = 4/3 pi r3
dShell[subSphere,'volume'] <- 4/3 * pi * (dShell[subSphere, 'xDim']/2 * dShell[subSphere, 'yDim']/2 * dShell[subSphere, 'zDim']/2)

#(mass^1/3)/volume^1/3 with shell volume substituted for shell size 
#do you need to take both to 1/3 sqrt?
dShell$densityMV <- dShell$mass1 / dShell$volume
#dShell[subWax,'volume'] <- NA
#dShell[subWax,'densityMV'] <- NA
#deviation from spherical and rootMass for an added bonus
pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3


#1.6 export this data as a csv
write.csv(dShell, 'shellsData_Expt2.csv', row.names=TRUE)

