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


#Script for analysis of cleaned data after processing
#Originally part of script 'Expt1 Analysis2'

setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt2.csv')
source('./PDF Prefs General.R')

  pData <- dShell
  TAXA <- unique(pData$taxon)
  pData$taxon <- as.factor(pData$taxon)
  taxa <- sort(unique(pData$taxon))
  tFont <- rep(3,length(taxa))
  tFont[which(taxa == 'Abranda')] <- 1
  taxaAbrev <- substring(taxa,0,5)

#0.0 Measurement checks general
par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(3,3,1,1))

  for (t in TAXA) {
    
    plot(xDim ~ yDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','yDim'))
    plot(xDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','zDim'))
    plot(yDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'yDim','zDim'))
    plot(xDim ~ thick, data=pData[(pData$taxon == t),], main=paste(t,'xDim','thick'))
    plot(xDim ~ mass1, data=pData[(pData$taxon == t),], main=paste(t,'xDim','mass1'))
    
    
  }

dev.off()

#2 Dissolution by total mass (mg, cMass) relationship to 4 variables: surface area (mm^2), initial mass, volume, and densityMV
#need to separate taxon out by color rather than letter

pdf('./outFigs/cMassMeasures.pdf', page='A4', height = 6 , width = 8)
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(4,4,1,1), cex=1)

  #2.0 first, mass lost by taxon 
  plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  
  #2.1 Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
  plot(cMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
  fit <- lm(cMass~finalSA, data=dShell)
  co <- coef(fit)
  abline(fit, lwd=2)
  plot(fit)
  summary(fit)
  
  #2.2 second, plot the mass lost/SA by taxon
  if (length(which(!is.na(pData$finalSA))) > 0) {
    plot(cMass/finalSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
    points(cMass/finalSA ~ taxon, data=pData)
    mtext('Mass lost / Surface Area (mg/mm\u00b2)', side=2, line=3)
    axis(2, las=1)
    axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
    #   mtext(e)
  } else {
    plot(1:1, type='n', ann=FALSE, axes=FALSE)
    mtext('no data', side=1, line=-2)
  }
  
  #2.3 plot mass by initial mass
  plot(cMass~mass1, data=pData, axes=TRUE, ann=TRUE)
  fit <- glm(cMass~mass1, data=dShell)
  co <- coef(fit)
  abline(fit, lwd=2)
  
  #2.4 plot volume
  plot(cMass~volume, data=pData, axes=TRUE, ann=TRUE)
  fit <- glm(cMass~volume, data=dShell)
  co <- coef(fit)
  abline(fit, lwd=2)
  
  plot(sqrt(cMass)~sqrt(volume), data=pData)
  fit <- glm(sqrt(cMass)~sqrt(volume), data=dShell)
  co <- coef(fit)
  abline(fit, lwd=2)
  
  #2.5 plot density
  plot(cMass~densityMV, data=pData, axes=TRUE, ann=TRUE)
  fit <- glm(cMass~densityMV, data=dShell)
  co <- coef(fit)
  abline(fit, lwd=2)

  #2.6 Mass lost/density by taxon
  if (length(which(!is.na(pData$densityMV))) > 0) {
    plot(cMass/densityMV ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
    points(cMass/densityMV ~ taxon, data=pData)
    mtext('Mass lost / density', side=2, line=3)
    axis(2, las=1)
    axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  } else {
    plot(1:1, type='n', ann=FALSE, axes=FALSE)
    mtext('no data', side=1, line=-2)
  }

dev.off()

#2.7 Modeling with CMass, continued

  #1 cMass with finalSA - decent linear relationship (when transformed)
  a5 <- plot(cMass ~ finalSA, data=pData)
  a5 <- lm(cMass ~ finalSA, data=pData)
  plot(a5)
  
  b5 <- plot(sqrt(cMass) ~ sqrt(finalSA), data=pData)
  b5 <- lm(sqrt(cMass) ~ sqrt(finalSA), data=pData)
  plot(b5)
  
  #2 cMass with mass1 - decent linear relationship (again when transformed)
  a6 <- plot(cMass ~ mass1, data=pData)
  a6 <- lm(cMass ~ mass1, data=pData)
  plot(a6)
  
  b6 <- plot(log(cMass) ~ log(mass1), data=pData)
  b6 <- lm(log(cMass) ~ log(mass1), data=pData)
  plot(b6)
  
  #3 cMass with volume - less decent linear relationship but still okay (when transformed)
  a7 <- plot(cMass ~ volume, data=pData)
  a7 <- lm(cMass ~ volume, data=pData)
  plot(a7)
  
  b7 <- plot(log(cMass) ~ log(volume), data=pData)
  b7 <- lm(log(cMass) ~ log(volume), data=pData)
  plot(b7)

  
#3.1 pMass: general plots for visualization

  #pdf setup
  pdf('./outFigs/pMassMeasures.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1))
  
    #1 percent mass lost by taxon
    plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='% Mass Lost')
    points(pMass ~ taxon, data=pData)
    mtext('% Mass Lost', side=2, line=3)
    axis(2, las=1)
    axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
    
    #2. percent mass lost by surface area
    plot(pMass ~ finalSA, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)')
    fit <- glm(pMass~finalSA, data=dShell)
    co <- coef(fit)
    abline(fit, lwd=2)
  
  dev.off() 
  
  #3 pMass with mass1 - somehow is an exponential relationship
  
  a2 <- plot(pMass ~ mass1, data=pData)
  a2 <- lm(pMass ~ mass1, data=pData)
  plot(a2)
  
  b2 <- plot(sqrt(pMass) ~ sqrt(mass1), data=pData)
  b2 <- lm(sqrt(pMass) ~ sqrt(mass1), data=pData)
  summary(b2)
  plot(b2)
  
  #4 pMass with densityMV - almost has potential, then fitted residuals end up being crap
  
  a4 <- plot(pMass ~ densityMV, data=pData)
  a4 <- lm(pMass ~ densityMV, data=pData)
  plot(a4)
  
  b4 <- plot(log(pMass) ~ log(densityMV), data=pData)
  b4 <- lm(log(pMass)~ log(densityMV), data=pData)
  plot(b4)
  