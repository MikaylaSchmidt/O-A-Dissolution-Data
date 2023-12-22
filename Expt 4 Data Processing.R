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


#Script for processing of raw data into full dataset used in analysis
#Originally part of script 'Expt1 Data Processing'


#0.open script, set names and parameters
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShellFull <- read.csv('./shells_Expt4.csv', skip=23)


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



#1.1 Surface area calculations for all of the taxa using caliper measurements
#by default, all are classified as cylinder
dShell$shape <- 'cylinder1'
dShell[(dShell$taxon == 'Liloa'),'shape'] <-'cylinder2'
dShell[(dShell$taxon == 'Calcite'),'shape'] <- 'cube'
dShell[(dShell$taxon == 'Tridacna'),'shape'] <- 'cube'

#turbo is dome + circle aka hemisphere
dShell[(dShell$taxon == 'Turbo'), 'shape'] <- 'hemisphere'

#pingui and abra technically two domes
dShell[(dShell$taxon == 'Pinguitellina'), 'shape'] <- '2dome'
dShell[(dShell$taxon == 'Abranda'), 'shape'] <- '2dome'
#while scaph is two cones
dShell[(dShell$taxon == 'Scaphopod'), 'shape'] <- '2cone'

#halimeda different shape entirely... diamond/rhombus
dShell[(dShell$taxon =='Halimeda'), 'shape'] <- 'rhombus'

#finally, nat and eth are spherical
dShell[(dShell$taxon == 'Natica'), 'shape'] <-'sphere'
dShell[(dShell$taxon == 'Ethalia'), 'shape'] <-'sphere'



#1.2 surface area caliper calculations
#basic rectangular cube surface area. does not work for halimeda based off how measurements were taken.
#is 2xy + 2yz + 2xz
dShell$cSA1 <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#next lines for surface area of cylindrical specimen, for margi cylinder
#equal to 2piR^2 + 2*h*pi*r
subCyl1 <- which(dShell$shape == 'cylinder1')
dShell[subCyl1,'cSA1'] <- 2 * (pi* dShell[subCyl1,'xDim']/2 * dShell[subCyl1,'yDim']/2) + dShell[subCyl1,'zDim']*pi*(dShell[subCyl1,'yDim'] + dShell[subCyl1,'xDim'])/2
#for liloa cylinder, 2piR^2 + 2*h*pi*r
subCyl2 <- which(dShell$shape == 'cylinder2')
dShell[subCyl2,'cSA1'] <- 2 * (pi* dShell[subCyl2,'yDim']/2 * dShell[subCyl2,'yDim']/2) + 2* dShell[subCyl2,'zDim']*pi*dShell[subCyl2,'yDim']/2 

#next turbo, a hemisphere
#is pi*r^2 + 2*pi*r*h + 
subHemi <- which(dShell$shape == 'hemisphere')
dShell[subHemi,'cSA1'] <-  pi * (dShell[subHemi,'xDim']/2 * dShell[subHemi,'yDim']/2) + 2 * pi * (dShell[subHemi,'xDim']/2 + dShell[subHemi,'yDim']/2)/2 * dShell[subHemi,'zDim'] 

#next lines for surface area of domed specimens
# is 2 * (2*pi*r*h)
subDome <- which(dShell$shape == '2dome')
dShell[subDome,'cSA1'] <- 2 * pi * 2 *(dShell[subDome,'xDim']/2 + dShell[subDome,'yDim']/2)/2 * dShell[subDome,'zDim']

#cone for scaph
# is 2 * pi*r(r+sqrt(h^2 +r^2)) 
subCone <- which(dShell$shape == '2cone')
dShell[subCone,'cSA1'] <- 2 * pi * dShell[subCone,'xDim']/2 *((dShell[subCone,'xDim']/2)+sqrt((dShell[subCone,'xDim']/2)^2+dShell[subCone,'zDim']^2))

#then, rough surface area for rhomboid Halimeda
#surface area is 2(xy/2) + 4(zc) where c = ((x/2)^2) + (y/2)^2)^1/2 
subRhom <- which(dShell$shape == 'rhombus')

#following calculation needs to be thrown out in favor of the imageJ calc
dShell[subRhom,'cDim'] <- sqrt((dShell[subRhom,'xDim']/2)^2 + (dShell[subRhom,'yDim']/2)^2)
dShell[subRhom,'cSA1'] <- 2 * ((dShell[subRhom,'xDim']*dShell[subRhom,'yDim'])/2) + 4 * dShell[subRhom,'zDim'] * dShell[subRhom,'cDim']


#finally, nat and eth spherical calc
#is 4 * pi * [(x/2 + y/2 + z/2)/3]^2
subSphere <- which(dShell$shape == 'sphere')
dShell[subSphere,'cSA1'] <- 4 * pi * ((dShell[subSphere,'xDim']/2+dShell[subSphere,'yDim']/2 + dShell[subSphere,'zDim']/2)/3)^2



#1.3 General calculation of surface area for each taxa based on ImageJ calculations
#Rough surface area for rhomboid Halimeda
#surface area is 2(xy/2) + 4(zc) where c = ((x/2)^2) + (y/2)^2)^1/2 
subRhom <- which(dShell$shape == 'rhombus')
dShell[subRhom,'calcSAPerim'] <- dShell[subRhom, 'perimeter'] * dShell[subRhom, 'zDim'] + dShell[subRhom, 'calcSA1'] + dShell[subRhom, 'calcSA2']


#quick comparison of default ImageJ SA (calcSA) with perimeter based SA (calcSA2) 
plot(calcSA~calcSAPerim, data=dShell)
abline(a= 0, b= 1)
test <- lm(calcSA~calcSAPerim, data=dShell)
summary(test)
#based off of linear model, it is significant to use perimeter calculated SA over approximation
dShell[subRhom,'calcSA'] <- dShell[subRhom,'calcSAPerim']

#1.5 Final dataset for surface area using alternate for halimeda. could also just use calcSA.
dShell$finalSA <- dShell$cSA1 
dShell[subRhom, 'finalSA'] <- dShell[subRhom, 'calcSAPerim'] 



#1.4 Density as given in Kosnik et al 2009

#Calculations of volume based on assigned shape
#cube volume = lwh 
dShell$volume <- dShell$xDim * dShell$yDim * dShell$zDim
#cylinder1 (margi) = pi r2 h  
dShell[subCyl1,'volume'] <- pi * (dShell[subCyl1, 'xDim']/2 * dShell[subCyl1, 'yDim']/2) * dShell[subCyl1, 'zDim']
#cylinder2 (liloa) = pi r2 h
dShell[subCyl2,'volume'] <- pi * (dShell[subCyl2, 'yDim']/2)^2 * dShell[subCyl2, 'zDim']
#dome (turbo) = 1/3 pi h^2 (3* r-h)
dShell[subHemi,'volume'] <- 1/3 * pi * dShell[subHemi,'zDim']^2 * (3 * (dShell[subHemi,'xDim']/2 + dShell[subHemi,'yDim']/2)/2 - dShell[subHemi,'zDim']) 
#2dome = dome1 full volume = 1/3 pi h^2 (3* r-h) 
dShell[subDome,'volume'] <- 1/3 * pi * dShell[subDome,'zDim']^2 * (3 * (dShell[subDome,'xDim']/2 + dShell[subDome,'yDim']/2)/2 - dShell[subDome,'zDim'])
#2cone = cone, volume = pi*r^2*(h/3)
dShell[subCone,'volume'] <- pi * (dShell[subCone,'xDim']/2)^2 * (dShell[subCone,'zDim']/3)
#rhombus(halimeda) = l w h with length * width being = the image J calculated surface area
dShell[subRhom,'volume'] <- (dShell[subRhom, 'calcSA1']+dShell[subRhom, 'calcSA2'])/2 * dShell[subRhom, 'zDim']
#sphere (nat and eth) = 4/3 pi r3
dShell[subSphere,'volume'] <- 4/3 * pi * (dShell[subSphere, 'xDim']/2 * dShell[subSphere, 'yDim']/2 * dShell[subSphere, 'zDim']/2)

#(mass^1/3)/volume^1/3 with shell volume substituted for shell size 
#do you need to take both to 1/3 sqrt?
dShell$densityMV <- dShell$cMass / dShell$volume



#1.6 export this data as a csv
write.csv(dShell, 'shellsData_Expt4.csv', row.names=TRUE)