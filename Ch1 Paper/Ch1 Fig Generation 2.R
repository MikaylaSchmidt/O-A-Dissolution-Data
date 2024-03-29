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
#pData <- subset(dShell, dShell$exptID == 'T1.1'|dShell$exptID == 'T1.2'|dShell$exptID == 'T1.3'|dShell$exptID == 'T4.1'|dShell$exptID == 'T4.2'|dShell$exptID == 'T4.3')
pData <- dShell
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

#few trial plots
plot(finalSA~cSize, data= pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2))
plot(rootMass~cSize, data= pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2))

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
legend('topright', legend= paste(TAXA2$taxon,' slope diff=',relSlope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#figure generation for shell cSize to mass1 ratio. should actually be done in terms of densityMV
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
plot(thick~cSize, data=tempNo, col= tempNo$tColor, pch=substring(tempNo$taxon, 0, 2), ylab='Thickness (mm)', xlab='Shell Size (mm)', ylim=c(0,12), xlim=c(0,12))
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

library(relevance)
TAXA3 <- dropdata(TAXA2, rowid= 'Calcite', incol='taxon')
TAXA3 <- dropdata(TAXA3, rowid='Tridacna', incol='taxon')
summaryData2 <- dropdata(summaryData, rowid= 'Calcite', incol='taxon.taxon')
summaryData2 <- dropdata(summaryData2, rowid='Tridacna', incol='taxon.taxon')

legend('topleft', legend= paste(TAXA3$taxon,' r\u00b2=',summaryData2$r2), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
legend('topleft', legend= paste(TAXA3$taxon,' slope =',summaryData2$slope), col=TAXA3$tColor, lty=TAXA3$tLine, text.font = TAXA3$tFont, pch= TAXA3$tPoint, text.col= TAXA3$tColor)
#remove calcite and tridacna from this because their thickness is arbitrary


#third one, for deviation from the 'normal' sphere
plot(deviation~cSize, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Deviation from Sphere (mm)', xlab='Shell Size (mm)', xlim=c(0,12))
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

#final one - for csize and SA by taxon
#not sure if needed or how to flip
plot(cSize ~ taxon, data=pData[!is.na(pData$cSize),], ylim=c(0,max(pData$cSize, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(cSize ~ taxon, data=pData[!is.na(pData$cSize),])
mtext('Mean Geometric Size (mm)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

plot(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),], ylim=c(0,max(pData$finalSA, na.rm=TRUE)), ann=FALSE, axes=FALSE)
points(finalSA ~ taxon, data=pData[!is.na(pData$finalSA),])
mtext('Surface Area', side=2, line=3)
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

#summary and residuals plotting, first OLS then SMA
summary(lm(percentTotal~percentSATotal, data = Expt123))
plot(lm(percentTotal~percentSATotal, data = Expt123))


#plot this
#but in pdf form
pdf('RatiosLost.pdf', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1))
plot(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2), xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentSATotal, data = Expt123)
abline(lm(percentTotal~percentSATotal, data = Expt123))
legend('topleft', legend= paste(TAXA2$taxon), col=TAXA2$tColor, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


abline(a=0, b=1)
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
print(summaryData)


#3.5 plot of total predictors
#plot with all of them included EXCEPT volume...
par(mfrow=c(1,1))
model0 <- lm(percentTotal ~ percentSATotal, data = Expt123)
summary(model0)
plot(model0)

Expt123$resid <- residuals(model0)
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of lm(percentTotal ~ percentSATotal)')
abline(h=0)

#issue with next is that densityMV and deviation are from same (?)
model1<- lm(percentTotal ~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo, data=Expt123)
summary(model1)
plot(model1)


model2<- lm(percentTotal ~ percentSATotal + cSize + thick + densityMV + deviation + exptID + exptNo + polymorph, data=Expt123)
summary(model2)
Expt123$resid <- residuals(model2)
plot(density(Expt123$resid))
boxplot(resid~taxon, data=Expt123, xlab = 'Taxon', ylab ='Residuals', main='Residuals of model2')
abline(h=0)

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
library(AICcmodavg)
AICc(model0, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model1, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model2, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model3, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model4, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model5, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model6, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model7, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model8, return.K = FALSE, second.ord = FALSE, nobs=NULL)
AICc(model9, return.K = FALSE, second.ord = FALSE, nobs=NULL)

useBIC(model0, return.K = FALSE, nobs=NULL)
useBIC(model1, return.K = FALSE, nobs=NULL)
useBIC(model2, return.K = FALSE, nobs=NULL)
useBIC(model3, return.K = FALSE, nobs=NULL)
useBIC(model4, return.K = FALSE, nobs=NULL)
useBIC(model5, return.K = FALSE, nobs=NULL)
useBIC(model6, return.K = FALSE, nobs=NULL)
useBIC(model7, return.K = FALSE, nobs=NULL)
useBIC(model8, return.K = FALSE, nobs=NULL)
useBIC(model9, return.K = FALSE, nobs=NULL)

#3.6 then do SMA analysis with these same plots to see if the results end up being better
library(smatr)
smaModel1 <- sma(percentTotal ~ mass1 + thick + finalSA + volume + densityMV, data=Expt123)
smaModel1

#issue with this specific one with dropped levels of grouping
smaModel2 <- sma(percentTotal ~ densityMV + finalSA, data=Expt123)
smaModel2 

smaModel3 <- sma(percentTotal ~ finalSA + densityMV + polymorph, data=Expt123)
smaModel3
smaModel4 <- sma(percentTotal~finalSA + densityMV + polymorph + taxon, data=Expt123)
smaModel4

lmodel2(percentTotal~percentSATotal, data = Expt123)
plot(lmodel2(percentTotal~percentSATotal, data = Expt123))
perTrial <- lmodel2(percentTotal~percentSATotal, data = Expt123) 
summary(perTrial)

perTrial2<-sma(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2))
plot(perTrial2, type = 'l', col= 'black', xlab = '% of Total Surface Area', ylab = "% of Total Mass Lost")
points(percentTotal~percentSATotal, data = Expt123, col= Expt123$tColor, pch=substring(Expt123$taxon, 0, 2))
