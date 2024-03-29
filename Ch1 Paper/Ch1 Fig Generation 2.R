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

#additional variables transformed
pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3


#figure generation for shell cSize to mass1 ratio
plot(rootMass~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
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
selectTemp <- which(pData$taxon %in% c('Tridacna','Calcite') == FALSE)
tempNo <- pData[(selectTemp),]
plot(thick~cSize, data=tempNo, col= tempNo$tColor, pch=substring(tempNo$taxon, 0, 2), ylab='Thickness', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
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
  summaryData[p,'n'] <- nrow(temp)
  if ((T %in% c('Tridacna','Calcite')) == FALSE)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
print(summaryData)

TAXA3 <- TAXA2[3:11,]
summaryData2 <- summaryData[3:11,]

legend('topleft', legend= paste(TAXA3$taxon,' r\u00b2=',summaryData2$r2), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
legend('topleft', legend= paste(TAXA3$taxon,' slope =',summaryData2$slope), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
#remove calcite and tridacna from this because their thickness is arbitrary

#third one, for deviation from the 'normal' sphere
plot(deviation~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Deviation', xlab='Shell Size', xlim=c(0,12))
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
abline(a=0, b=1, lty=2, lwd=1.55)
abline(a=0, b=0, lty=2, lwd=1.55)
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
legend('topleft', legend= paste(TAXA2$taxon,' slope =',1-slope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

print(summaryData)

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

pData$rootMass <- (pData$mass1)^ (1/3)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3

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
        
            #Work on this
            stat_boxplot()
            table <- A$stats
            

#just the two controls
#you could redo these looking at individual experiments rather than the average
subControl <- subset(Expt1_3, Expt1_3$taxon == 'Calcite'|Expt1_3$taxon == 'Tridacna')
ggplot(subControl, aes(taxon, avg)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic() + xlab('Taxon') + ylab('% Dissolution per Hour') +
  scale_fill_discrete(name = "Starting pH", labels = c("pH 5.1", "pH 7.1"))


dev.off()


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

#summary and residuals plotting, first OLS then SMA
summary(lm(percentTotal~percentSATotal, data = Expt1_3))
plot(lm(percentTotal~percentSATotal, data = Expt1_3))


#plot this
#but in pdf form
pdf('RatiosLost.pdf', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1))
plot(percentTotal~percentSATotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2), xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentSATotal, data = Expt1_3)
abline(lm(percentTotal~percentSATotal, data = Expt1_3))
legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#abline(a=0, b=1)
#summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
#summary(summaryData)
#p = 0
#for(T in TAXA2$taxon){
  #temp <- Expt1_3[(Expt1_3$taxon == T),]
  #pref <- TAXA2[(TAXA2$taxon==T),]
  #lm.temp <- lm(percentTotal~percentSATotal, data=temp)
  #p <- p + 1
  #summaryData[p,'r2'] <- round(summary(lm.temp)$adj.r.squared,3)
  #summaryData[p, 'p'] <- round(summary(lm.temp)$coefficients[2,4],7)
  #summaryData[p,'slope'] <- round(summary(lm.temp)$coefficients[2,1],3)
  #summaryData[p,'slope.err'] <- round(summary(lm.temp)$coefficients[2,2],3)
  #summaryData[p,'b.est'] <- round(summary(lm.temp)$coefficients[1,1],3)
  #summaryData[p,'b.err'] <- round(summary(lm.temp)$coefficients[1,2], 3)
  #summaryData[p,'n'] <- nrow(temp)
  #abline(lm.temp, col=pref$tColor, lty=pref$tLine)
#}


#3.5 plot of total predictors
#plot with all of them included EXCEPT volume...
par(mfrow=c(1,1))
model0 <- lm(percentTotal ~ percentSATotal, data = Expt1_3)
summary(model0)

#issue with next is that densityMV and deviation are from same (?)
model1<- lm(percentTotal ~ thick + percentSATotal + densityMV + deviation + rootMass + mass1 + exptNo + exptID, data=Expt1_3)
summary(model1)

model2 <- lm(percentTotal ~ thick + percentSATotal + densityMV + deviation + polymorph + rootMass, data=Expt1_3)
summary(model2)

model3 <- lm(percentTotal~ percentSATotal + densityMV + deviation + polymorph, data=Expt1_3)
summary(model3)

model4<- lm(percentTotal~percentSATotal + densityMV + deviation + polymorph, data=Expt1_3)
summary(model4)
plot(model4)

model5<- lm(percentTotal~ percentSATotal + thick  + densityMV + deviation  + polymorph + taxon, data=Expt1_3)
summary(model5)
plot(model5)



#3.6 then do SMA analysis with these same plots to see if the results end up being better
library(smatr)
smaModel1 <- sma(percentTotal ~ mass1 + thick + finalSA + volume + densityMV, data=Expt1_3)
smaModel1

#issue with this specific one with dropped levels of grouping
smaModel2 <- sma(percentTotal ~ densityMV + finalSA, data=Expt1_3)
smaModel2 

smaModel3 <- sma(percentTotal ~ finalSA + densityMV + polymorph, data=Expt1_3)
smaModel3
smaModel4 <- sma(percentTotal~finalSA + densityMV + polymorph + taxon, data=Expt1_3)
smaModel4

lmodel2(percentTotal~percentSATotal, data = Expt1_3)
plot(lmodel2(percentTotal~percentSATotal, data = Expt1_3))
perTrial <- lmodel2(percentTotal~percentSATotal, data = Expt1_3) 
summary(perTrial)

perTrial2<-sma(percentTotal~percentSATotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2))
plot(perTrial2, type = 'l', col= 'black', xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
points(percentTotal~percentSATotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2))
