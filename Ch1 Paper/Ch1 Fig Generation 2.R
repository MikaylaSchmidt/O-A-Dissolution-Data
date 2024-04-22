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
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

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

#additional variables transformed
pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3

#first, looking at csize and SA by taxon
#not sure if needed or how to flip
plot(cSize ~ taxon, data=pData[!is.na(pData$cSize),], ylim=c(0,max(pData$cSize, na.rm=TRUE)), ann=FALSE, axes= FALSE )
points(cSize ~ taxon, data=pData[!is.na(pData$cSize),])
mtext('Geometric Mean Size (mm)', side=2, line=3)
axis(2, las=1)
abline(h=4)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

plot(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),], ylim=c(0,max(pData$finalSA, na.rm=TRUE)), xlab= '', xaxt= 'n', ann=FALSE, yaxt='n')
points(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),])
mtext('Surface Area (mm\U00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
mtext('Taxon', side = 1, line = 6)

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
pData$poros <- pData$mass1/(pData$volume * 2.83)
subCal <- which(pData$polymorph =='Calcite')
pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
pData$poros <- 1- pData$poros
plot(poros~taxon, data=pData, ylab = 'Porosity (Actual Mass/Expected Mass)', xlab ='', xaxt='n' )
points(poros~taxon, data = pData)
abline(h=0)
mtext('Taxon', side = 1, line = 6)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

#figure generation for shell cSize to mass1 ratio. should actually be done in terms of densityMV
plot(rootMass~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Specimen Size', ylim=c(2,11), xlim=c(2,15))
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

pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3

#3.1 Combining Expt 1, Expt 2, and Expt 3 and then subsetting by experiment
#Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
pData$poros <- pData$mass1/(pData$volume * 2.83)
subCal <- which(pData$polymorph =='Calcite')
pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
Expt123 <- pData
Expt123$exptNo <- 'Expt1'
Expt123[(Expt123$exptID == 'T4.1'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.2'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T4.3'),'exptNo'] <-'Expt3'
Expt123[(Expt123$exptID == 'T2.1'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.2'),'exptNo'] <-'Expt2'
Expt123[(Expt123$exptID == 'T2.3'),'exptNo'] <-'Expt2'

factor(Expt123$exptNo)
Expt123$exptNo <- as.factor(Expt123$exptNo)
str(Expt123)

#total carbonate per replicate in mg
cTotal <- aggregate(Expt123$cMass, by = list(Expt123$exptID),FUN=sum)
colnames(cTotal)<-c('exptID','cTotal')
Expt123 <- merge(Expt123, cTotal, by = 'exptID')

#3.2 for plotting
#% calcium carbonate dissolved over total
Expt123$percentTotal <- (Expt123$cMass)/(Expt123$cTotal) *100

#now for surface area
SATotal <- aggregate(Expt123$finalSA, by = list(Expt123$exptID),FUN=sum)
colnames(SATotal)<-c('exptID','SATotal')
Expt123 <- merge(Expt123, SATotal, by = 'exptID')
Expt123$percentSATotal <- (Expt123$finalSA)/(Expt123$SATotal) * 100


# NEED TO MOVE THIS BELOW EXPT 123 SUBSET
#anova test (parametric) and kruskal test (nonparametric)
#for taxon by surface area two different ways
res.aov <- aov(percentSATotal ~ taxon, data = Expt123)
summary(res.aov)
res.aov <- aov(pMass ~ taxon, data = Expt123)
summary(res.aov)
kruskal.test(percentSATotal ~ taxon, data= Expt123)
kruskal.test(pMass ~ taxon, data= Expt123)

#either % SA Total or final SA works - trends by taxon still the same
plot(percentSATotal ~ taxon, data=Expt123[!is.na(Expt123$percentSATotal),], ylim=c(0,max(Expt123$percentSATotal, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(percentSATotal ~ taxon, data=Expt123[!is.na(Expt123$percentSATotal),])
mtext('Surface Area (mm\U00b2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

res.aov <- aov(percentSATotal ~ taxon, data = Expt123)
summary(res.aov)


#trial plot of sqrtfinalSA by cSize
ExptNoWax <- subset(Expt123, waxYN == 'No Wax')
plot(sqrt(finalSA)~cSize, data= ExptNoWax, col= pData$tColor, pch=substring(pData$taxon, 0, 2))
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- ExptNoWax[(ExptNoWax$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(sqrt(finalSA)~cSize, data=temp)
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
legend('topright', legend= paste(TAXA2$taxon,' slope =',summaryData$slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)



#summary and residuals plotting, first OLS then SMA
summary(lm(percentTotal~percentSATotal, data = Expt123))
par(mfrow=c(2,2))
plot(lm(percentTotal~percentSATotal, data = Expt123))

#not exactly normally distributed 
#but that's not something you're going to get when you have data start at 0 
par(mfrow=c(1,1))
library(ggplot2)
library(ggpubr)
ggdensity(Expt123$percentTotal)
shapiro.test(Expt123$percentTotal)

#plot this
#but in pdf form
#pdf('RatiosLost.pdf', height = 6 , width = 8)
#par(mfrow=c(1,1), oma=c(1,1,1,1))
#Expt123 <- subset(Expt123, Expt123$taxon !='Calcite')
plot(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
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
  lm.temp <- lm(percentTotal~percentSATotal, data=temp)
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

#plot of what the residuals are demonstrating - some have more or less
par(mfrow=c(1,1))
Expt123$perRatio <- Expt123$percentTotal/Expt123$percentSATotal
plot(perRatio~taxon, data=Expt123)
abline(h=1)

#3.5 plot of total predictors
#plot with all of them included EXCEPT volume...
modelDef <- lm(percentTotal ~ percentSATotal, data = Expt123)
summary(modelDef)
plot(modelDef)
dredge(modelDef)

Expt123$resid <- residuals(modelDef)
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of lm(percentTotal ~ percentSATotal)')
abline(h=0)

#first, no taxon version for model selection
#create new data frame using only the relevant data aka relData
relData <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'poros', 'deviation', 'exptID', 'exptNo', 'polymorph','taxon')]
noTaxon <- relData[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'poros', 'deviation', 'exptID', 'exptNo', 'polymorph')]

#dredge function
library(lmerTest)
library(MuMIn)
model0 <- lm(percentTotal ~., data= relData)
options(na.action = 'na.fail')
summary(model0)
AICc(model0)
dredge(model0, evaluate = TRUE)
modelBest0 <- lm(percentTotal ~ percentSATotal + thick + cSize + taxon, data= relData)
dredge(modelBest0)
summary(modelBest0) 

model1 <- lm(percentTotal ~., data= noTaxon)
dredge(model1)
summary(model1)
AICc(model1)
modelBest1 <- lm(percentTotal ~ percentSATotal + cSize + polymorph, data= noTaxon)
summary(modelBest1)
AICc(modelBest1)

#linear mixed effects model doesn't change major takeaways
#need to use conditional aic instead for model selection
library(cAIC4)
model0.1 <- lmer(percentTotal~percentSATotal + cSize + thick + poros + deviation + polymorph + (1|exptID) + (1|taxon), data = relData)
dredge(model0.1)
summary(model0.1)
cAIC(model0.1)
r.squaredGLMM(lmer(percentTotal~percentSATotal + cSize + thick + poros + deviation + polymorph + (1|exptID) + (1|taxon), data = relData))
modelBest0.1 <- lmer(percentTotal~percentSATotal + polymorph + (1|taxon), data = relData)
rand(modelBest0.1)
r.squaredGLMM(lmer(percentTotal~percentSATotal + polymorph + (1|taxon), data = relData))
cAIC(modelBest0.1)
AIC(modelBest0.1)
cAIC(model0)

#next bit is honestly just a repeat of the last thing
model0.2 <- lmer(Expt123$percentTotal ~ Expt123$percentSATotal + Expt123$poros + Expt123$thick + Expt123$deviation + Expt123$cSize + Expt123$polymorph + (1|Expt123$exptID) + (1|Expt123$taxon))
dredge(model0.2)
step(model0.2)
rand(model0.2)
r.squaredGLMM(lmer(percentTotal~percentSATotal + cSize + polymorph + (1|taxon), data = relData))


#modelSelect <- lm(relData$percentTotal ~ relData$percentSATotal + relData$cSize + relData$densityMV + relData$thick + relData$deviation + relData$polymorph + relData$taxon + relData$exptID + relData$exptNo)
#linear mixed effects model
summary(lmer(log(Expt123$cMass)~ log(Expt123$finalSA) + log(Expt123$mass1) + log(Expt123$volume) + log(Expt123$thick) + Expt123$deviation + log(Expt123$cSize) + Expt123$polymorph + (1|Expt123$taxon) + (1|Expt123$exptID)))
rand(lmer(log(Expt123$cMass)~ log(Expt123$finalSA) + log(Expt123$mass1) + log(Expt123$volume) + log(Expt123$thick) + Expt123$deviation + log(Expt123$cSize) + Expt123$polymorph + (1|Expt123$taxon) + (1|Expt123$exptID)))
r.squaredGLMM(lmer(log(Expt123$cMass)~ log(Expt123$finalSA) + log(Expt123$mass1) + log(Expt123$volume) + log(Expt123$thick) + Expt123$deviation + log(Expt123$cSize) + Expt123$polymorph + (1|Expt123$exptID) + (1|Expt123$exptNo)))

#mixed effects with %
summary(lmer(Expt123$percentTotal~ Expt123$percentSATotal + Expt123$densityMV + Expt123$thick + Expt123$deviation + Expt123$cSize + Expt123$polymorph + (1|Expt123$exptID)))
rand(lmer(Expt123$percentTotal~ Expt123$percentSATotal + Expt123$densityMV + Expt123$thick + Expt123$deviation + Expt123$cSize + Expt123$polymorph + (1|Expt123$taxon) + (1|Expt123$exptID) + (1|Expt123$exptNo)))
r.squaredGLMM(lmer(Expt123$percentTotal~ Expt123$percentSATotal + Expt123$densityMV + Expt123$thick + Expt123$deviation + Expt123$cSize + Expt123$polymorph + (1|Expt123$taxon) + (1|Expt123$exptID) + (1|Expt123$exptNo)))


#issue with next is that densityMV and deviation are from same (?)
modelFull <- lm(percentTotal ~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo + polymorph + taxon, data=Expt123)
summary(modelFull)

model1<- lm(percentTotal ~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo, data=Expt123)
summary(model1)
plot(model1)


model2<- lm(percentTotal ~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo + polymorph, data=Expt123)
summary(model2)
Expt123$resid <- residuals(model2)
plot(density(Expt123$resid))
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of model2')
abline(h=0)
AIC(model2)

model3 <- lm(percentTotal~ percentSATotal + cSize + deviation + polymorph, data=Expt123)
summary(model3)

model4<- lm(percentTotal~percentSATotal + densityMV + deviation + polymorph + cSize, data=Expt123)
summary(model4)
plot(model4)

model5<- lm(percentTotal~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo + taxon, data=Expt123)
summary(model5)

model6 <- lm(percentTotal~ percentSATotal + cSize + thick + densityMV + deviation + polymorph + taxon, data=Expt123)
summary(model6)
plot(model6)

model7 <- lm(percentTotal~ percentSATotal + cSize + polymorph + taxon, data=Expt123)
summary(model7)

model8 <- lm(percentTotal~ percentSATotal + taxon, data=Expt123)
summary(model8)

model9 <- lm(percentTotal~ percentSATotal + polymorph, data=Expt123)
summary(model9)

#time for some AIC comparisons

#logit function lm
#new relData
library(boot)
logitData <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'thick', 'densityMV', 'deviation', 'exptID', 'exptNo', 'polymorph','taxon')]
logitData$percentTotal <- logit(logitData$percentTotal/100)
logitData$percentSATotal <- logit(logitData$percentSATotal/100)
logitModel <- lm(percentTotal ~., data = logitData)
plot(logitModel)
step(logitModel, direction ='backward')

#lm of logit transform surface area/mass lost (check diagnostic plots)
logitModel2 <- lm(percentTotal ~ percentSATotal, data = logitData)
plot(logitModel2)
