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


#Script for analysis of cleaned data after processing, for Ch2 Expt 1
#0.0 Setup for everything
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell <- read.csv('./Ch2Clean_Expt1.csv')

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Anadara')] <- 1
taxaAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')


#1.0 First, cMass relationships  
  #Which is mass lost by taxon 
  plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  
  #Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
  plot(cMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
  fit <- lm(cMass~finalSA, data=pData)
  co <- coef(fit)
  abline(fit, lwd=2)
  plot(fit)
  summary(fit)


#2.0 Then, looking at pMass relationships
  #percent lost by taxon 
  plot(pMass ~ taxon, data=pData[!is.na(pData$pMass),], ylim=c(0,max(pData$pMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(pMass ~ taxon, data=pData[!is.na(pData$pMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  
  #Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
  plot(pMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (%)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
  fit <- lm(pMass~finalSA, data=pData)
  co <- coef(fit)
  abline(fit, lwd=2)
  plot(fit)
  summary(fit)
  
#3.0 two plots looking at the deviation from spherical
  #third one, for deviation from the 'normal' sphere
  plot(deviation~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Deviation from Sphere (mm)', xlab='Shell Size (mm)', xlim=c(2,12))
  summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
  summary(summaryData)
  p = 0
  for(T in TAXA2$taxon){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(deviation~cSize, data=temp)
    p <- p + 1
    summaryData[p,'r2'] <- round(summary(lm.temp)$adj.r.squared,3)
    summaryData[p, 'p'] <- round(summary(lm.temp)$coefficients[2,4],7)
    summaryData[p,'slope'] <- round(summary(lm.temp)$coefficients[2,1],3)
    summaryData[p,'slope.err'] <- round(summary(lm.temp)$coefficients[2,2],3)
    summaryData[p,'b.est'] <- round(summary(lm.temp)$coefficients[1,1],3)
    summaryData[p,'b.err'] <- round(summary(lm.temp)$coefficients[1,2], 3)
    summaryData[p,'n'] <- nrow(temp)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
  }
  abline(a=0, b=1, lty=2, lwd=1.60)
  abline(a=0, b=0, lty=2, lwd=1.60)
  legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',summaryData$r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  legend('topleft', legend= paste(TAXA2$taxon,' slope diff=',1-summaryData$slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  
  print(summaryData)
  
  #alternative deviation plot
  plot(deviation~taxon, data=pData, xlab ='', xaxt ='n', ylab = 'Devation from Spherical (mm)')
  points(deviation~taxon, data = pData)
  abline(h=0)
  mtext('Taxon', side = 1, line = 6)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  
  
#4.0 densityMV, which should be substituted for the csize metric because is more accurate metric
  plot(mass1~volume, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Starting Mass (mg)', xlab='Volume (mm\U00b3)')
  arag <- data.frame(volume=0:1500,mass1=(0:1500)*2.83)
  lines(arag$volume, arag$mass1)
  lm(mass1 ~ volume, data=arag)
  calcite <- data.frame(volume=0:1500,mass1=(0:1500)*2.711)
  lines(calcite$volume, calcite$mass1, lty='dashed')
  lm(volume ~ mass1, data=calcite)
  summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
  summary(summaryData)
  p = 0
  for(T in TAXA2$taxon){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(mass1~volume, data=temp)
    p <- p + 1
    summaryData[p,'r2'] <- round(summary(lm.temp)$adj.r.squared,3)
    summaryData[p, 'p'] <- round(summary(lm.temp)$coefficients[2,4],7)
    summaryData[p,'slope'] <- round(summary(lm.temp)$coefficients[2,1],3)
    summaryData[p,'slope.err'] <- round(summary(lm.temp)$coefficients[2,2],3)
    summaryData[p,'b.est'] <- round(summary(lm.temp)$coefficients[1,1],3)
    summaryData[p,'b.err'] <- round(summary(lm.temp)$coefficients[1,2], 3)
    summaryData[p,'n'] <- nrow(temp)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
  }

  
#5.0 now getting into the actual comparison bit - looking at MassSA
  plot(massSA ~ taxon, data=pData)
  points(massSA ~ taxon, data=pData)
  mtext('Mass lost/surface area (mg/mm2)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2) 
  
  #Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
  #Same as seen in this script previously but is related to above massSA calc
  plot(cMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
  fit <- lm(cMass~finalSA, data=pData)
  co <- coef(fit)
  abline(fit, lwd=2)
  plot(fit)
  summary(fit)
 
#6.0 changing it into %/% standard
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
  plot(perMassSA~taxon, data= pData)
  points(perMassSA~taxon, data = pData)  
  
  #plots: how much is actually good
  plot(percentTotal~percentSATotal, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2))
  try<-lm(percentTotal~percentSATotal, data=pData)
  summary(try)
  try2<-lm(cMass ~ finalSA, data=pData)
  summary(try2)
  
  
  
#7.0 Comparing stuff - what if you look at just calcite specs? Does that improve things???????
  #first, calcite subset
  onlyCalc <- subset(pData, pData$polymorph == 'Calcite')
  droplevels(as.factor(onlyCalc$polymorph))
  droplevels(as.factor(onlyCalc$taxon))
  onlyCalc$taxon <- factor(onlyCalc$taxon, levels = c('Calcite',  'Centrostephanus', 'Pecten', 'Saccostrea'))
  
  
  #plot of massSA
  boxplot(massSA ~ taxon, data = onlyCalc)
  points(massSA~taxon, data=onlyCalc)
  
  #is mass/SA plot better fit by only calcite?
  plot(cMass ~ finalSA, data=onlyCalc)
  only1 <- lm(cMass ~ finalSA, data=onlyCalc)
  summary(only1)
  only2 <- lm(cMass ~ finalSA * mag, data=onlyCalc)
  summary(only2)
  
  #plot of perMassSA
  boxplot(perMassSA ~taxon, data = onlyCalc)
  points(perMassSA~taxon, data=onlyCalc)
       
  #perMassSA plot fit with just calcite 
  plot(percentTotal ~ percentSATotal, data=onlyCalc)
  only3 <- lm(percentTotal ~ percentSATotal, data=onlyCalc)
  summary(only3)
  only4 <- lm(percentTotal ~ percentSATotal * mag, data=onlyCalc)
  summary(only4)
  
#8.0 Comparing stuff using t tests
#asking the question: are these populations significantly different than each other?
  
  