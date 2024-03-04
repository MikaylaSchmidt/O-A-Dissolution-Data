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

#1.0 Script for creation of figures for expt 1
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')

#setup
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')
head(pData)


#1.1 graph 1 pMass
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 Paper Figs")
pdf('./Expt1pMass.pdf', width=8, height=6)
par(mfrow=c(1,1), oma=c(1,1,1,1))

#combined from script 'Expt1 Post Processing Analysis'
  #weirdly exponential 
  plot(mass1~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Starting Mass (mg)', xlab='% Mass Lost')
  pVal <- vector()
  p = 0
  for(T in TAXA){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(mass1~pMass, data=temp)
    p <- p + 1
    pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
  }
  legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  
  #less of a trend with % compared to surface area,try square root in order to linearity
  plot(finalSA~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Surface Area (mm\u00b2)', xlab='% Mass Lost')
  pVal <- vector()
  p = 0
  for(T in TAXA){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(finalSA~pMass, data=temp)
    p <- p + 1
    pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
  }
  legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  
  #more of a trend with % compared to density
  plot(densityMV~pMass, data=pData, col= pData$tColor , pch=substring(pData$taxon, 0, 2), ylab='Density (mg/mm\u00b3)', xlab='% Mass Lost')
  lm.den <- lm(densityMV~pMass, data=pData)
  summary(lm.den)
  lm.den1 <- lm(densityMV~pMass+taxon, data=pData)
  pVal <- vector()
  p = 0
  summary(lm.den1)
  for(T in TAXA){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(densityMV~pMass, data=temp)
    p <- p + 1
    pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
    }
  legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  #legend('topleft', legend= TAXA2$taxon, col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  
  #data with thickness is simply weird, but seems to related to density
  #need to remove cal/tridacna data because thickness arbitrary
  plot(thick~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Thickness (mm)', xlab='% Mass Lost')
  pVal <- vector()
  p = 0
  for(T in TAXA){
    temp <- pData[(pData$taxon == T),]
    pref <- TAXA2[(TAXA2$taxon==T),]
    lm.temp <- lm(thick~pMass, data=temp)
    p <- p + 1
    pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
    abline(lm.temp, col=pref$tColor, lty=pref$tLine)
  }
  legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
  
  
  #plot with all of them included EXCEPT volume...
  par(mfrow=c(1,1), oma=c(1,1,1,1))
  total <- pData$mass1 + pData$finalSA + pData$densityMV + pData$thick + pData$volume
  plot(total~pMass, data=pData, pch=substring(pData$taxon, 0, 2), ylab='Variables', xlab='% Mass Lost')
  
dev.off()


#1.2 graphs of total mass and % mass lost by taxon
pdf('./Expt1byTaxon.pdf', width=8, height=6)
par(mfrow=c(1,1), oma=c(1,1,1,1))

  #first with pMass
  plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='% Mass Lost')
  points(pMass ~ taxon, data=pData)
  mtext('% Mass Lost', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #then with cMass 
  plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  
dev.off()
  


#2.0 repeated, just for Expt 3 (formerly Expt 4)
#setup script
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt4.csv')

#more setup
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')
head(pData)

#2.1 graph 1 pMass
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 Paper Figs")
pdf('./Expt3pMass.pdf', width=8, height=6)
par(mfrow=c(1,1), oma=c(1,1,1,1))

#combined from script 'Expt1 Post Processing Analysis'
#weirdly exponential 
plot(mass1~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Starting Mass (mg)', xlab='% Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(mass1~pMass, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#less of a trend with % compared to surface area
plot(finalSA~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Surface Area (mm\u00b2)', xlab='% Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(finalSA~pMass, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#more of a trend with % compared to density
plot(densityMV~pMass, data=pData, col= pData$tColor , pch=substring(pData$taxon, 0, 2), ylab='Density (mg/mm\u00b3)', xlab='% Mass Lost')
lm.den <- lm(densityMV~pMass, data=pData)
summary(lm.den)
lm.den1 <- lm(densityMV~pMass+taxon, data=pData)
pVal <- vector()
p = 0
summary(lm.den1)
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(densityMV~pMass, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
#legend('topleft', legend= TAXA2$taxon, col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#data with thickness is simply weird, but seems to related to density
plot(thick~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Thickness (mm)', xlab='% Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(thick~pMass, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#plot with all of them included EXCEPT volume...
par(mfrow=c(1,1), oma=c(1,1,1,1))
total <- pData$mass1 + pData$finalSA + pData$densityMV + pData$thick + pData$volume
plot(total~pMass, data=pData, pch=substring(pData$taxon, 0, 2), ylab='Variables', xlab='% Mass Lost')

dev.off()

#2.2 comparisons by taxon
pdf('./Expt3byTaxon.pdf', width=8, height=6)
par(mfrow=c(1,1), oma=c(1,1,1,1))

  #first with pMass
  plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='% Mass Lost')
  points(pMass ~ taxon, data=pData)
  mtext('% Mass Lost', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #then with cMass 
  plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

dev.off()


