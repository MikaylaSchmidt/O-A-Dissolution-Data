################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2023 Mikayla
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

dShell <- read.csv('./shellsData_Expt1.csv')

#1.1 Comparison of the Caliper vs ImageJ measurements
#Setup for measurement comp  by taxa
pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Aragonite')] <- 1
taxaAbrev <- substring(taxa,0,5)

#1.2 Plot of caliper SA / ImageJ SA with abline at 1
plot(cSA1/calcSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='', xlim=c(0,11))
points(cSA1/calcSA ~ taxon, data=pData)
abline(a=1, b=0)
mtext('Caliper SA / Calculated SA', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxaAbrev, font=tFont, cex=0.5, las=2)

#1.3 Log plot of SAs
pData2 <- pData[!is.na(pData$calcSA),]
plot(cSA1~calcSA, data=pData2, pch=substring(pData2$taxon, 0, 2), log='xy', xlab='log(ImageJSA)', ylab='log(CaliperSA)')
abline(a=0, b=1)
levels(pData$taxon)

#1.4 Log plot of SAs by taxon
plot((cSA1-calcSA)/calcSA~taxon, data=pData2, pch=substring(pData$taxon, 0, 2))
log2(((pData2$cSA1-pData2$calcSA)/pData2$calcSA)+1)

#1.5 Mann-Whitney U Tests of each taxon (for ImageJ to caliper measures)
#Subset all 11 taxa for tests

Abra <- subset(dShell, taxon =='Abranda')
Abra1<-wilcox.test(Abra$cSA1, Abra$calcSA, data=dShell)
print(Abra1)

Arag <- subset(dShell, taxon =='Aragonite')
Arag1<-wilcox.test(Arag$cSA1, Arag$calcSA, data=dShell)
print(Arag1)

Eth <- subset(dShell, taxon =='Ethalia')
Eth1<-wilcox.test(Eth$cSA1, Eth$calcSA, data=dShell)
print(Eth1)

Hali <- subset(dShell, taxon =='Halimeda')
Hali1<-wilcox.test(Hali$cSA1, Hali$calcSA, data=dShell)
print(Hali1)

Liloa <- subset(dShell, taxon =='Liloa')
Liloa1<-wilcox.test(Liloa$cSA1, Liloa$calcSA, data=dShell)
print(Liloa1)

Margi <- subset(dShell, taxon =='Marginopora')
Margi1<-wilcox.test(Margi$cSA1, Margi$calcSA, data=dShell)
print(Margi1)

Nat <- subset(dShell, taxon =='Natica')
Nat1<-wilcox.test(Nat$cSA1, Nat$calcSA, data=dShell)
print(Nat1)

Pingui <- subset(dShell, taxon =='Pinguitellina')
Pingui1<-wilcox.test(Pingui$cSA1, Pingui$calcSA, data=dShell)
print(Pingui1)

Scaph <- subset(dShell, taxon =='Scaphopod')
Scaph1<-wilcox.test(Scaph$cSA1, Scaph$calcSA, data=dShell)
print(Scaph1)

Tri <- subset(dShell, taxon =='Tridacna')
Tri1<-wilcox.test(Tri$cSA1, Tri$calcSA, data=dShell)
print(Tri1)

Turbo <- subset(dShell, taxon =='Turbo')
Turbo1<-wilcox.test(Turbo$cSA1, Turbo$calcSA, data=dShell)
print(Turbo1)


#2.1 Dissolution (mg) relationship to surface area (mm^2)
#need to separate taxon out by color rather than letter

plot(cMass ~ cSA1, data=dShell, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (mg)',pch=substring(dShell$taxon, 0, 2))
fit <- glm(cMass~cSA1, data=dShell)
co <- coef(fit)
abline(fit, lwd=2)


#2.2 second, plot the mass vs calc by taxon

plot(cMass/calcSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
points(cMass/calcSA ~ taxon, data=pData)
mtext('Mass lost / Calculated SA', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)



#2.3 Modeling this all - start with basic linear models
#Do you need to change pMass (%) to cMass (total in mg) 
a <- lm(pMass ~ calcSA, data=pData)
summary(a)
b <- lm(pMass ~ calcSA + mass1, data=pData)
summary(b)
c <- lm(pMass ~ calcSA + taxon, data=pData)
summary(c)
d <- lm(pMass ~ calcSA + mass1 + taxon, data=pData)
summary(d)
e <- lm(pMass ~ calcSA + density, data=pData)
summary(e)
f <- lm(pMass ~ calcSA + mass1 + density, data=pData)
summary(f)
g <- lm(pMass ~ calcSA + mass1 + taxon + density, data=pData)
summary(g)

a$residuals

#don't forget to add model looking at replicates




#3.1 Comparison of the pH over time 
#read in data
pH1.1 <- read.csv('./Expt 1.1 Log 1.csv')
pH1.2 <- read.csv('./Expt 1.2 Log 1.csv')
pH1.3 <- read.csv('./Expt 1.3 Log 1.csv')

#convert interval to hours for each (interval # x 5 min / 60)
pH1.1$hours <- (pH1.1$Interval * 5)/60 
pH1.2$hours <- (pH1.2$Interval * 5)/60 
pH1.3$hours <- (pH1.3$Interval * 5)/60 

#time ~ pH for replicate 1.1
plot(pH ~ hours, data=pH1.1, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.1 pH Over Time")
fit1.1 <- glm(pH ~ hours, data=pH1.1)
co <- coef(fit1.1)
abline(fit1.1, lwd=2)

#repeat for 1.2
plot(pH ~ hours, data=pH1.2, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.2 pH Over Time")
fit1.2 <- glm(pH ~ hours, data=pH1.2)
co <- coef(fit1.2)
abline(fit1.2, lwd=2)

#and 1.3
plot(pH ~ hours, data=pH1.3, type = "o", xlab = "Time (hours)", ylab = "pH", main = "Replicate 1.3 pH Over Time")
fit1.3 <- glm(pH ~ hours, data=pH1.3)
co <- coef(fit1.3)
abline(fit1.3, lwd=2)

dev.off()

#3.2 Examining pH 
pH1 <- lm(pH ~ hours, data=pH1.1)
summary(pH1)
pH2 <- lm(pH ~ hours, data=pH1.2)
summary(pH2)
pH3 <- lm(pH ~ hours, data=pH1.3)
summary(pH3)
