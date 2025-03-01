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
dShellFull <- read.csv('./Ch2_Expt1.csv', skip=22)


#remove missing or broken taxon, remove two bad columns
dShell <- subset(dShellFull, exclude == 0)
dShell <- dShell[-c(14:15)]
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
  
  #this works for oyster specimens
  #is 2xy + 2yz + 2xz
  dShell$finalSA <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

  #updated surface area for pecten - two wiggly sides
  subPect <- which(dShell$taxon == 'Pecten')
  dShell[subPect,'calcX'] <- (dShell[subPect,'calcX1'] + dShell[subPect,'calcX2'])/2
  dShell[subPect,'finalSA'] <- 2 * dShell[subPect,'calcX'] * dShell[subPect,'yDim'] + 2 * dShell[subPect,'yDim'] * dShell[subPect,'zDim'] + 2 * dShell[subPect,'calcX'] * dShell[subPect,'zDim']
  
  #updated surface area for anadara - one wiggly side
  subAna <- which(dShell$taxon == 'Anadara')
  dShell[subAna,'calcX'] <- (dShell[subAna,'calcX1'] + dShell[subAna,'calcX2'])/2
  dShell[subAna,'finalSA'] <- dShell[subAna,'calcX'] * dShell[subAna,'yDim'] + dShell[subAna,'xDim'] * dShell[subAna,'yDim'] + 2 * dShell[subAna,'yDim'] * dShell[subAna,'zDim'] + 2 * dShell[subAna,'calcX'] * dShell[subAna,'zDim']
  
  
  #And for cone...
  #only using the equation for the lateral portion of the cone, as the base contains another cone
  #is 2 * (pi * r * l) with l being slanted side of the cone
  #aka 2 * pi*r(sqrt(h^2 +r^2)) 
  subCone <- which(dShell$taxon == 'Centrostephanus')
  dShell[subCone,'finalSA'] <- 2 * pi *dShell[subCone,'zDim']/2*(sqrt((dShell[subCone,'zDim']/2)^2+dShell[subCone,'xDim']^2))
  dShell$massSA <- dShell$cMass/dShell$finalSA
  


#2.0 volume and density calculations
#Calculations of volume based on assigned shape

  #cube volume = lwh 
  dShell$volume <- dShell$xDim * dShell$yDim * dShell$zDim
  
  #cube volume for the two complicated cubes
  dShell[subPect,'volume'] <- dShell[subPect,'calcX'] * dShell[subPect,'yDim'] * dShell[subPect,'zDim']
  dShell[subAna, 'volume'] <- dShell[subAna, 'xDim'] * dShell[subAna,'yDim'] * dShell[subAna,'zDim']
  
  #for cone, volume = pi*r^2*(h/3)
  dShell[subCone,'volume'] <- pi * (dShell[subCone,'yDim']/2 * dShell[subCone,'zDim']/2) * (dShell[subCone,'xDim']/3)
  
  #(mass^1/3)/volume^1/3 with shell volume substituted for shell size 
  dShell$densityMV <- dShell$mass1 / dShell$volume
  
  #deviation from spherical and rootMass for an added bonus
  dShell$rootMass <- (dShell$mass1)^ (1/3)
  dShell$deviation <- (abs(dShell$xDim - dShell$cSize) + abs(dShell$yDim - dShell$cSize) +  abs(dShell$zDim - dShell$cSize))/3
  dShell[subPect, 'deviation'] <- (abs(dShell[subPect,'calcX'] - dShell[subPect,'cSize']) + abs(dShell[subPect,'yDim'] - dShell[subPect,'cSize']) +  abs(dShell[subPect,'zDim'] - dShell[subPect,'cSize']))/3
  dShell$deviationStan <- dShell$deviation/dShell$cSize
  

#3.0 calculation of porosity using calcite 
  dShell$poros <- dShell$mass1/(dShell$volume * 2.711)
  subCal <- which(dShell$polymorph =='Aragonite')
  dShell[subCal, 'poros'] <- dShell[subCal, 'mass1'] / (dShell[subCal, 'volume'] * 2.93)
  dShell[subAna, 'poros'] <- dShell[subAna, 'mass1'] / (dShell[subAna, 'volume'] * 2.8862)
  dShell$poros <- 1- dShell$poros



#4.0 export this data as a csv
write.csv(dShell, 'Ch2Clean_Expt1.csv', row.names=TRUE)
dev.off()

