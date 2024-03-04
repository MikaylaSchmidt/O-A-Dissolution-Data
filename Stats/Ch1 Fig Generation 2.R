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

#read in data first
#1.0 Script for creation of figures for ch 1
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

#setup in order to prep data & exclude Trial 2 experimental data (wax makes things bad)
pData <- subset(dShell, dShell$exptID == 'T1.1'|dShell$exptID == 'T1.2'|dShell$exptID == 'T1.3'|dShell$exptID == 'T4.1'|dShell$exptID == 'T4.2'|dShell$exptID == 'T4.3')
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

#figure generation for shell cSize to mass1 ratio
pData$rootMass <- (pData$mass1)^ (1/3)
plot(rootMass~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
arag <- data.frame(cSize=0:12,rootMass=(0:12)*2.83)
arag$rootMass = ((arag$cSize^3)*2.83)^(1/3)
lines(arag$cSize, arag$rootMass)
lm(rootMass ~ cSize, data=arag)
calcite <- data.frame(cSize=0:12,rootMass=(0:12)*2.711)
calcite$rootMass = ((calcite$cSize^3)*2.711)^(1/3)
lines(calcite$cSize, calcite$rootMass)
lm(rootMass ~ cSize, data=calcite)
r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(rootMass~cSize, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}

#getting p values
#work on this later
pVal <- vector()
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(rootMass~cSize, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$p.value,3)

}
summary(lm.temp)
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#comment legend below is original r^2 legend, while beyond that is for slope
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

relSlope <- round(slope/1.41,3)
relSlope
cbind(TAXA2,relSlope)
legend('topleft', legend= paste(TAXA2$taxon,' slope=',relSlope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#same but for thickness
selectTemp <- which(pData$taxon %in% c('Tridacna','Calcite') == FALSE)
tempNo <- pData[(selectTemp),]
plot(thick~cSize, data=tempNo, col= tempNo$tColor, pch=substring(tempNo$taxon, 0, 2), ylab='Thickness', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(thick~cSize, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  if ((T %in% c('Tridacna','Calcite')) == FALSE)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
#remove calcite and tridacna from this because their thickness is arbitrary

#third one, for deviation from the 'normal' sphere
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3
plot(deviation~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Deviation', xlab='Shell Size', xlim=c(0,12))
r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(deviation~cSize, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
abline(h=1)
abline(a=0, b=1)
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',1-slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#final one - for csize by taxon
#not sure if needed or how to flip
plot(cSize ~ taxon, data=pData[!is.na(pData$cSize),], ylim=c(0,max(pData$cSize, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(cSize ~ taxon, data=pData[!is.na(pData$cSize),])
mtext('Mean Geometric Size', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)


#2.0 now a comparison of the two, dependent on measurement
#setup script
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')


#3.1 Combining Expt 1 and Expt 3 without Expt 2, and then subsetting by experiment

Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
Expt1_3$exptNo <- 'Expt1'
Expt1_3[(Expt1_3$exptID == 'T4.1'),'exptNo'] <-'Expt3'
Expt1_3[(Expt1_3$exptID == 'T4.2'),'exptNo'] <-'Expt3'
Expt1_3[(Expt1_3$exptID == 'T4.3'),'exptNo'] <-'Expt3'


#3.2 Setup and execution of comparison plots
#only thing necessary for this is to split it into two more subsets by overall experiment 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 Paper Figs")
pdf('Expt1&3Comp.pdf', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1))

#firstly, some setup to subset
library(ggplot2)
factor(Expt1_3$exptNo)
Expt1_3$exptNo <- as.factor(Expt1_3$exptNo)
str(Expt1_3)

#cMass by taxon, split by experiment
ggplot(Expt1_3, aes(taxon, cMass)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic()

#pMass by taxon, split by experiment
ggplot(Expt1_3, aes(taxon, pMass)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic()

#pMass/hour by taxon, split by experiment
#first, general avg
Expt1_3$avg <- ((Expt1_3$pMass)/48) * 100

#then for expt3
sub3 <- which(Expt1_3$exptNo == 'Expt3')
Expt1_3[sub3,'avg'] <- (Expt1_3[sub3,'pMass']/168) *100
ggplot(Expt1_3, aes(taxon, avg)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic() + xlab('Taxon') + ylab('% Dissolution per Hour') +
  scale_fill_discrete(name = "Starting pH", labels = c("pH 5.1", "pH 7.1"))

#just the two controls
#you could redo these looking at individual experiments rather than the average
subControl <- subset(Expt1_3, Expt1_3$taxon == 'Calcite'|Expt1_3$taxon == 'Tridacna')
ggplot(subControl, aes(taxon, avg)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic() + xlab('Taxon') + ylab('% Dissolution per Hour') +
  scale_fill_discrete(name = "Starting pH", labels = c("pH 5.1", "pH 7.1"))


dev.off()



#3.2.2 shell mean dissolution by taxon
Expt1_3$meanDiss <- ((Expt1_3$cMass)^(2/3))/48 
sub3 <- which(Expt1_3$exptNo == 'Expt3')
Expt1_3[sub3,'meanDiss'] <- Expt1_3[sub3,'cMass']^(2/3)/168

plot(meanDiss ~ taxon, data=Expt1_3, ann=FALSE, axes=FALSE)
points(meanDiss ~ taxon, data=Expt1_3)
mtext('Mean Dissolution per Hour', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

plot(meanDiss~mass1, data=Expt1_3, col= Expt1_3$tColor , pch=substring(Expt1_3$taxon, 0, 2))
r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA){
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(meanDiss~mass1, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)




#3.3 comparison plots but standardized for pH
#subsets to add ph in to dataset, calculated separately.
#H+ is measured in molarity (mol/liter)
Expt1_3$deltaH <- 0.0000060473
Expt1_3[(Expt1_3$exptID == 'T1.2'),'deltaH'] <-0.0000042753
Expt1_3[(Expt1_3$exptID == 'T1.3'),'deltaH'] <-0.0000043809
Expt1_3[(Expt1_3$exptID == 'T4.1'),'deltaH'] <-0.0000000262
Expt1_3[(Expt1_3$exptID == 'T4.2'),'deltaH'] <-0.0000000488
Expt1_3[(Expt1_3$exptID == 'T4.3'),'deltaH'] <-0.0000000832

#change mol/liter hydrogen to mol CaCO3
Expt1_3$predictCarb <- ((Expt1_3$deltaH * 30.1)/2)
sub3 <- which(Expt1_3$exptNo == 'Expt3')
Expt1_3[sub3,'predictCarb'] <- (Expt1_3[sub3,'deltaH'] *30.075)/2

#total carbonate per replicate in mg
cTotal <- aggregate(Expt1_3$cMass, by = list(Expt1_3$exptID),FUN=sum)
colnames(cTotal)<-c('exptID','cTotal')
Expt1_3 <- merge(Expt1_3, cTotal, by = 'exptID')

#cMass is currently in mg, needs to be in mol like H+
#below equation goes from calcium carbonate solid mg to g to mol (divided by liters)
Expt1_3$cTotalMol <- (Expt1_3$cTotal)/(100.0869*1000)

#plot predicted carbonate over actual carbonate in mols
plot(predictCarb~cTotalMol, data= Expt1_3)

#3.4 for plotting
#% calcium carbonate dissolved over total
Expt1_3$percentTotal <- (Expt1_3$cMass)/(Expt1_3$cTotal) *100

#now for surface area
SATotal <- aggregate(Expt1_3$finalSA, by = list(Expt1_3$exptID),FUN=sum)
colnames(SATotal)<-c('exptID','SATotal')
Expt1_3 <- merge(Expt1_3, SATotal, by = 'exptID')
Expt1_3$percentSATotal <- (Expt1_3$finalSA)/(Expt1_3$SATotal) * 100

summary(lm(percentTotal~percentSATotal, data = Expt1_3))

#plot this
#but in pdf form
pdf('RatiosLost.pdf', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1))
plot(percentTotal~percentSATotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2), xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentSATotal, data = Expt1_3)
abline(lm(percentTotal~percentSATotal, data = Expt1_3))
abline(a=0, b=1)

r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(percentTotal~percentSATotal, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#again but for volume instead of SA 
volumeTotal <- aggregate(Expt1_3$volume, by = list(Expt1_3$exptID),FUN=sum)
colnames(volumeTotal)<-c('exptID','volumeTotal')
Expt1_3 <- merge(Expt1_3, volumeTotal, by = 'exptID')
Expt1_3$percentVolumeTotal <- (Expt1_3$volume)/(Expt1_3$volumeTotal) * 100

plot(percentTotal~percentVolumeTotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2), xlab = '% of Total Volume', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentVolumeTotal, data = Expt1_3)
abline(lm(percentTotal~percentVolumeTotal, data = Expt1_3))
abline(a=0, b=1)

r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(percentTotal~percentVolumeTotal, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

dev.off()

#SA seems to be 1:1 while volume is closer to 1:2
Expt1_3$newRatio <- (Expt1_3$percentTotal)/(Expt1_3$percentVolumeTotal)
plot(newRatio~taxon, data=Expt1_3, ylab='% Mass Lost/% Volume')
abline (a=1, b=0)
