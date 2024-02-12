#read in data first
#1.0 Script for creation of figures for expt 1
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
calcite$rootMass = ((arag$cSize^3)*2.711)^(1/3)
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
#comment legenth below is original r^2 legend, while beyond that is for slope
#legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
relSlope <- round(slope/1.41,3)
relSlope
cbind(TAXA2,relSlope)
legend('topleft', legend= paste(TAXA2$taxon,' slope=',relSlope), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#same but for thickness
#using TAXA3 data, which doesn't include calcite or tridacna/aragonite
selectTemp <- which(pData$taxon %in% c('Tridacna','Calcite') == FALSE)
tempNo <- pData[(selectTemp),]

plot(thick~cSize, data=pData, col= tempNo$tColor, pch=substring(tempNo$taxon, 0, 2), ylab='Thickness', xlab='Shell Size', ylim=c(0,12), xlim=c(0,12))
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
summary(lm.temp)

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
abline(h=0)
abline(a=1, b=0)
legend('topleft', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
