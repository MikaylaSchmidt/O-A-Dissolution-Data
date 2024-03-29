#FMA/SMA analysis

#default OLS lm here, using just R
test1<-lm(rootMass~cSize, data=pData)
summary(test1)

#comparative package for OLS/SMA/RMA
install.packages("lmodel2")

#read in lmodel2 package
library(lmodel2)
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

#first bit for rootmass by csize
pData$rootMass <- (pData$mass1)^ (1/3)
lmodel2(rootMass~cSize, data=pData)
rootTrial <- lmodel2(rootMass~cSize, data=pData)
plot(rootTrial, 'SMA') 
plot(rootTrial, 'SMA', col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))

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

#trying to do the same thing again but with sma for each individual taxon
rootTrial <- sma(rootMass~cSize, data=pData)
plot(rootTrial, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
summaryData <- data.frame(taxon=TAXA2, n=NA,  p=NA, r2=NA, b.est=NA, b.err=NA, slope=NA, slope.err=NA)
summary(summaryData)
p = 0
for(T in TAXA2$taxon){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  sma.temp <- sma(rootMass~cSize, data=temp)
  p <- p + 1
  summaryData[p,'r2'] <- round(summary(sma.temp)$adj.r.squared,3)
  summaryData[p, 'p'] <- round(summary(sma.temp)$coefficients[2,4],7)
  summaryData[p,'slope'] <- round(summary(sma.temp)$coefficients[2,1],3)
  summaryData[p,'slope.err'] <- round(summary(sma.temp)$coefficients[2,2],3)
  summaryData[p,'b.est'] <- round(summary(sma.temp)$coefficients[1,1],3)
  summaryData[p,'b.err'] <- round(summary(sma.temp)$coefficients[1,2], 3)
  summaryData[p,'n'] <- nrow(temp)
  abline(sma.temp, col=pref$tColor, lty=pref$tLine)
}


#or with alternative SMAtr package
install.packages("smatr")
library(smatr)
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

pData <- subset(dShell, dShell$exptID == 'T1.1'|dShell$exptID == 'T1.2'|dShell$exptID == 'T1.3'|dShell$exptID == 'T4.1'|dShell$exptID == 'T4.2'|dShell$exptID == 'T4.3')
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Calcite')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')

#first bit for rootmass by csize
pData$rootMass <- (pData$mass1)^ (1/3)
sma(rootMass~cSize, data=pData, method = 'SMA')
rootTrial <- sma(rootMass~cSize, data=pData)
plot(rootTrial, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Cube Root Mass', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))


