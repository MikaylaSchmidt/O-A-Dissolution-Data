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

#read in data first
#1.0 Script for creation of figures for ch 1
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./Ch1Clean_MasterSet.csv')
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

#setup in order to prep data & exclude Trial 2 experimental data (wax makes things bad)
#pData <- subset(dShell, dShell$exptID == 'T1.1'|dShell$exptID == 'T1.2'|dShell$exptID == 'T1.3'|dShell$exptID == 'T4.1'|dShell$exptID == 'T4.2'|dShell$exptID == 'T4.3')
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
TAXA2$taxon <- factor(TAXA2$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))
TAXA2 <- TAXA2[order(levels(TAXA2$taxon)),]
pData<-merge(pData,TAXA2,by='taxon')
pData$taxon <- factor(pData$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))
taxa <- sort(unique(pData$taxon))

#natica fix
#subNatNew <- which(pData$taxon == 'Notocochlis' & pData$waxYN == 'No Wax')
#pData[subNatNew,'finalSA']<- pData[subNatNew, 'finalSA'] * 1.5


#additional variables transformed
pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))
pData$deviation <- pData$deviation/pData$cSize

#first, looking at csize and SA by taxon
#not sure if needed or how to flip
plot(cSize ~ taxon, data=pData[!is.na(pData$cSize),], ylim=c(0,max(pData$cSize, na.rm=TRUE)), ann=FALSE, axes= FALSE )
points(cSize ~ taxon, data=pData[!is.na(pData$cSize),])
mtext('Geometric Mean Size (mm)', side=2, line=3)
axis(2, las=1)
abline(h=4)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#first this in native r
plot(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),], ylim=c(0,max(pData$finalSA, na.rm=TRUE)), xlab= '', xaxt= 'n', ann=FALSE, yaxt='n')
points(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),])
mtext('Surface Area (mm\U00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
mtext('Taxon', side = 1, line = 6)

#then in ggplot
ggplot(pData, aes(taxon, finalSA),) +
  geom_boxplot() + theme_clean() + labs(x = 'Taxon') +
  labs(y='Surface Area (mm\U00b2)')

plot(finalSA ~ cSize, data=pData[!is.na(pData$finalSA),], ylim=c(0,max(pData$finalSA, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(finalSA ~ cSize, data=pData[!is.na(pData$finalSA),])
mtext('Surface Area (mm\U00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#densityMV, which should be substituted for the csize metric because is more accurate metric
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
print(summaryData)

#comment legend below is original r^2 legend, while beyond that is for slope
#legend('topright', legend= paste(TAXA2$taxon,'r\u00b2=', summaryData$r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

relSlope <- round(summaryData$slope/2.83,3)
relSlope
cbind(TAXA2,relSlope)
legend('bottomright', legend= paste(TAXA2$taxon,' rel slope=',relSlope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#alternative porosity plot
#pData$poros <- pData$mass1/(pData$volume * 2.83)
#subCal <- which(pData$polymorph =='Calcite')
#pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
#pData$poros <- 1- pData$poros
plot(poros~taxon, data=pData, ylab = 'Porosity (Actual Mass/Expected Mass)', xlab ='', xaxt='n' )
points(poros~taxon, data = pData)
abline(h=0)
mtext('Taxon', side = 1, line = 6)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#figure generation for shell cSize to mass1 ratio. should actually be done in terms of densityMV
plot(rootMass~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Specimen Size')
arag <- data.frame(cSize=0:12,rootMass=(0:12)*2.83)
arag$rootMass = ((arag$cSize^3)*2.83)^(1/3)
lines(arag$cSize, arag$rootMass)
lm(rootMass ~ cSize, data=arag)
calcite <- data.frame(cSize=0:12,rootMass=(0:12)*2.711)
calcite$rootMass = ((calcite$cSize^3)*2.711)^(1/3)
lines(calcite$cSize, calcite$rootMass)
lm(rootMass ~ cSize, data=calcite)
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(rootMass~cSize, data=temp)
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

#update code to pull r2/slope from above 
summary(lm.temp)
print(summaryData)

#comment legend below is original r^2 legend, while beyond that is for slope
legend('topright', legend= paste(TAXA2$taxon,'r\u00b2=', summaryData$r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

relSlope <- round(summaryData$slope/1.41,3)
relSlope
cbind(TAXA2,relSlope)
legend('topleft', legend= paste(TAXA2$taxon,' slope=',relSlope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#same but for thickness
selectTemp <- which(pData$taxon %in% c('Aragonite','Calcite') == FALSE)
tempNo <- pData[(selectTemp),]
plot(thick~cSize, data=tempNo, col= tempNo$tColor, pch=substring(tempNo$taxon, 0, 2), ylab='Thickness (mm)', xlab='Specimen Size (mm)', ylim=c(0,5), xlim=c(2,12))
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA, df = NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(thick~cSize, data=temp)
  p <- p + 1
  summaryData[p,'r2'] <- round(summary(lm.temp)$adj.r.squared,3)
  summaryData[p, 'p'] <- round(summary(lm.temp)$coefficients[2,4],7)
  summaryData[p,'slope'] <- round(summary(lm.temp)$coefficients[2,1],3)
  summaryData[p,'slope.err'] <- round(summary(lm.temp)$coefficients[2,2],3)
  summaryData[p,'b.est'] <- round(summary(lm.temp)$coefficients[1,1],3)
  summaryData[p,'b.err'] <- round(summary(lm.temp)$coefficients[1,2], 3)
  summaryData[p,'df'] <- lm.temp$df.residual
  summaryData[p,'n'] <- nrow(temp)
  if ((T %in% c('Aragonite','Calcite')) == FALSE)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
abline(h=0, lwd=2)
print(summaryData)

library(relevance)
TAXA3 <- dropdata(TAXA2, rowid= 'Calcite', incol='taxon')
TAXA3 <- dropdata(TAXA3, rowid='Aragonite', incol='taxon')
summaryData2 <- dropdata(summaryData, rowid= 'Calcite', incol='taxon.taxon')
summaryData2 <- dropdata(summaryData2, rowid='Aragonite', incol='taxon.taxon')

#legend('topleft', legend= paste(TAXA3$taxon,' r\u00b2=',summaryData2$r2), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
droplevels(TAXA3$taxon)
droplevels(summaryData2$taxon.taxon)
legend('topright', legend= paste(TAXA3$taxon,' slope =',summaryData2$slope), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
#remove calcite and tridacna from this because their thickness is arbitrary


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



#2.0 now a comparison of the two, dependent on measurement
#setup script
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)
dShell <- read.csv('./Ch1Clean_MasterSet.csv')

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

#natica fix
#subNatNew <- which(pData$taxon == 'Notocochlis' & pData$waxYN == 'No Wax')
#pData[subNatNew,'finalSA']<- pData[subNatNew, 'finalSA'] * 1.5

#pData$rootMass <- (pData$mass1)^ (1/3)
#pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))
#pData$deviation <- pData$deviation/pData$cSize
#pData$poros <- pData$mass1/(pData$volume * 2.83)
#subCal <- which(pData$polymorph =='Calcite')
#pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)



#3.1 Combining Expt 1, Expt 2, and Expt 3 and then subsetting by experiment
#Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
Expt123 <- pData
Expt123$exptNo <- 'Expt1'
Expt123[(Expt123$exptID == 'T4.1'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.2'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.3'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T2.1'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.2'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.3'),'exptNo'] <-'Expt2'
#Expt123$massSA <- Expt123$cMass/Expt123$finalSA

factor(Expt123$exptNo)
Expt123$exptNo <- as.factor(Expt123$exptNo)
str(Expt123)

  #start of no calcite breakout, only run this if calculating aragonite only 
  Expt123 <- subset(Expt123, Expt123$taxon !='Calcite')
  Expt123 <- subset(Expt123, Expt123$taxon !='Marginopora')
  #Expt123 <- subset(Expt123, Expt123$taxon !='Turbo')
  Expt123$taxon <- factor(Expt123$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Aragonite'))
  

#total carbonate per replicate in mg
#cTotal <- aggregate(Expt123$cMass, by = list(Expt123$exptID),FUN=sum)
#colnames(cTotal)<-c('exptID','cTotal')
#Expt123 <- merge(Expt123, cTotal, by = 'exptID')

#3.2 for plotting
#% calcium carbonate dissolved over total
#Expt123$percentTotal <- (Expt123$cMass)/(Expt123$cTotal) *100

#now for surface area
#SATotal <- aggregate(Expt123$finalSA, by = list(Expt123$exptID),FUN=sum)
#colnames(SATotal)<-c('exptID','SATotal')
#Expt123 <- merge(Expt123, SATotal, by = 'exptID')
#Expt123$percentSATotal <- (Expt123$finalSA)/(Expt123$SATotal) * 100


#either % SA Total or final SA works - trends by taxon still the same
plot(percentSATotal ~ taxon, data=Expt123[!is.na(Expt123$percentSATotal),], ylim=c(0,max(Expt123$percentSATotal, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(percentSATotal ~ taxon, data=Expt123[!is.na(Expt123$percentSATotal),])
mtext('Standardized Surface Area (mm\U00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

res.aov <- aov(percentSATotal ~ taxon, data = Expt123)
summary(res.aov)

#loading in packages and looking at density/ normality
par(mfrow=c(1,1))
library(ggplot2)
library(ggpubr)
ggdensity(Expt123$percentTotal/Expt123$percentSATotal)
shapiro.test(Expt123$percentTotal/Expt123$percentSATotal)

#not exactly normally distributed 
#but that's not something you're going to get when you have data start at 0 
#now for the linear model itself
summary(lm(percentTotal~percentSATotal, data = Expt123))
par(mfrow=c(2,2))
plot(lm(percentTotal~percentSATotal, data = Expt123))

#redo taxa2 to not include margi/arag/turbo
TAXA2 <- subset(TAXA2, TAXA2$taxon !='Calcite')
TAXA2 <- subset(TAXA2, TAXA2$taxon !='Marginopora')
#TAXA2 <- subset(TAXA2, TAXA2$taxon !='Turbo')
TAXA2$taxon <- factor(TAXA2$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Aragonite'))

#quick plot of just cMass by SA
plot(cMass~finalSA, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = 'Surface Area', ylab = "Mass Lost")
linePlot <-lm(cMass~finalSA, data = Expt123)
summary(linePlot)
abline(lm(cMass~finalSA, data = Expt123))

#this should correlate with avg arag mass lost/SA
mean(Expt123$massSA)

#IMPORTANT LINEAR MODEL 
#plot this
#but in pdf form
#pdf('RatiosLost.pdf', height = 6 , width = 8)
par(mfrow=c(1,1))
plot(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = 'Standardized Experimental Surface Area', ylab = "Standardized Experimental Mass Lost")
line <-lm(percentTotal~percentSATotal, data = Expt123)
abline(lm(percentTotal~percentSATotal, data = Expt123))
legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
abline(a=0, b=1, lty='dashed')


summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- Expt123[(Expt123$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(percentTotal ~ percentSATotal, data=temp)
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
legend('topleft', legend= paste(TAXA2$taxon,'slope =',summaryData$slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

#summaryData$oneDiff <- 'Y'
#summaryData$lmDiff <- 'N'
print(summaryData)
summaryData$taxon.taxon <- factor(summaryData$taxon.taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Aragonite'))


#plot of what the residuals are demonstrating - some have more or less
par(mfrow=c(1,1))
Expt123$perRatio <- Expt123$percentTotal/Expt123$percentSATotal
plot(perRatio~taxon, data=Expt123)
abline(h=1)

#3.5 plot of total predictors
#plot with just the main percent experimental surface area...
modelDef <- lm(percentTotal ~ percentSATotal, data = Expt123)
summary(modelDef)
plot(modelDef)
dredge(modelDef)
AICc(modelDef)

Expt123$resid <- residuals(modelDef)
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of lm(percentTotal ~ percentSATotal)')
abline(h=0)

#above + taxon
modelDef2 <- lm(percentTotal ~ percentSATotal + taxon, data = Expt123)
summary(modelDef2)
plot(modelDef2)
dredge(modelDef2)
AICc(modelDef2)


#first, no taxon version for model selection
#create new data frame using only the relevant data aka relData
relData <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'poros', 'deviationStan', 'exptID','exptNo', 'taxon')]
noTaxon <- relData[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'poros', 'deviationStan', 'exptID' , 'exptNo')]

#correlation check for variables
corDataNum <- Expt123[,c('cMass', 'finalSA', 'cSize', 'thick', 'poros', 'deviation')]
cor(corDataNum)
relDataNum <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'poros', 'deviation')]
cor(relDataNum)

#dredge function
library(lmerTest)
library(MuMIn)
model1 <- lm(percentTotal ~., data= noTaxon)
options(na.action = 'na.fail')
dredge(model1)
summary(model1)
AICc(model1)
modelBest1 <- lm(percentTotal ~ percentSATotal + deviationStan, data= noTaxon)
summary(modelBest1)
AICc(modelBest1)

model0 <- lm(percentTotal ~., data= relData)
summary(model0)
AICc(model0)
dredge(model0, evaluate = TRUE)
modelBest0 <- lm(percentTotal ~ percentSATotal  + cSize + taxon, data= relData)
dredge(modelBest0)
summary(modelBest0) 
AICc(modelBest0)


#model with just size and then size + taxon
modelSize <- lm(percentTotal ~ cSize, data = Expt123)
summary(modelSize)
modelSizeTax <- lm(percentTotal ~ cSize + taxon, data = Expt123)
summary(modelSizeTax)

#linear mixed effects model doesn't change major takeaways
#need to use conditional aic instead for model selection
library(cAIC4)
model0.1 <- lmer(percentTotal~percentSATotal + cSize + thick + poros + deviation  + (1|exptID)+ (1|exptNo) + (1|taxon), data = relData)
dredge(model0.1)
summary(model0.1)
cAIC(model0.1)
r.squaredGLMM(lmer(percentTotal~percentSATotal + cSize + thick + poros + deviation + (1|taxon), data = relData))
modelBest0.1 <- lmer(percentTotal~percentSATotal + cSize + (1|taxon), data = relData)
rand(modelBest0.1)
r.squaredGLMM(lmer(percentTotal~percentSATotal + (1|taxon), data = relData))
cAIC(modelBest0.1)
AIC(modelBest0.1)
cAIC(model0)

#next bit is honestly just a repeat of the last thing
model0.2 <- lmer(Expt123$percentTotal ~ Expt123$percentSATotal + Expt123$poros + Expt123$thick + Expt123$deviation + Expt123$cSize + Expt123$polymorph + (1|Expt123$exptID) + (1|Expt123$taxon))
dredge(model0.2)
step(model0.2)
rand(model0.2)
r.squaredGLMM(lmer(percentTotal~percentSATotal + cSize + polymorph + (1|taxon), data = relData))


##additional calcite trial to have a look at those values
#don't run the calcite breakout for this one 
subCal <- subset(Expt123, Expt123$taxon== "Calcite")
plot(percentTotal~percentSATotal, data = subCal)
abline(lm(percentTotal~percentSATotal, data = subCal))

subCal1 <- subset(Expt123, Expt123$taxon== "Calcite" & Expt123$exptNo == "Expt1" )
subCal2 <- subset(Expt123, Expt123$taxon== "Calcite" & Expt123$exptNo == "Expt2" )
subCal3 <- subset(Expt123, Expt123$taxon== "Calcite" & Expt123$exptNo == "Expt3" )
plot(percentTotal~percentSATotal, data = subCal1, col= 'blue', xlim = c(0,4), ylim = c(0,3.5), xlab ='Percent of Experimental Surface Area', ylab = 'Percent of Experimental Mass Lost')
points(percentTotal ~percentSATotal, data = subCal2, col='red')
points(percentTotal ~ percentSATotal, data = subCal3, col = 'green')
abline(lm(percentTotal~percentSATotal, data = subCal1), col = 'blue')
abline(lm(percentTotal~percentSATotal, data = subCal2), col='red')
abline(lm(percentTotal~percentSATotal, data = subCal3), col='green')
leg.text <- c('Expt 1', 'Expt 2', 'Expt 3')
legend('topleft', legend=paste(leg.text), col = c('blue','red','green'), pch=1)

plot(cMass~finalSA, data = subCal1, col= 'blue', xlim = c(0,400), ylim = c(0,15), xlab ='Surface Area', ylab = 'Mass Lost')
points(cMass~finalSA, data = subCal2, col='red')
points(cMass~finalSA, data = subCal3, col = 'green')
abline(lm(cMass~finalSA, data = subCal1), col = 'blue')
abline(lm(cMass~finalSA, data = subCal2), col='red')
abline(lm(cMass~finalSA, data = subCal3), col='green')
leg.text <- c('Expt 1', 'Expt 2', 'Expt 3')
legend('topleft', legend=paste(leg.text), col = c('blue','red','green'), pch=1)

cal3Try <- lm(massSA ~ exptID, data=subCal3)
anova(cal3Try)

