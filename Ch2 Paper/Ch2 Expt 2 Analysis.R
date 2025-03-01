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


#Script for analysis of cleaned data after processing, for Ch2 Expt 2
#0.0 Setup for everything
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell <- read.csv('./Ch2Clean_Expt2.csv')

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Anadara')] <- 1
taxaAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
#TAXA2$taxon <- factor(TAXA2$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))
#TAXA2 <- TAXA2[order(levels(TAXA2$taxon)),]
pData<-merge(pData,TAXA2,by='taxon')
#pData$taxon <- factor(pData$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))
pData$massSA <- pData$cMass/pData$finalSA


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

#percent lost by taxon by sealer presence
pData$seal <- factor(pData$seal)
noDef <- subset(pData, pData$taxon != 'Aragonite')
noDef <- subset(noDef, noDef$taxon != 'Calcite')
noDef <- subset(noDef, noDef$taxon != 'Centrostephanus')
droplevels(noDef$taxon)
boxplot(pMass ~ seal + droplevels(noDef$taxon), data=noDef)
stripchart(pMass ~ seal + droplevels(noDef$taxon), data = noDef, vertical=TRUE, add=TRUE, pch=1)


#Dissolution by total mass (mg, cMass) relationship to surface area (mm^2)
plot(pMass ~ finalSA, data=pData, xlab='Surface Area (mm\u00b2)', ylab='Dissolution (%)', col= pData$tColor, pch=substring(pData$taxon, 0, 2))
fit <- lm(pMass~finalSA, data=pData)
co <- coef(fit)
abline(fit, lwd=2)
plot(fit)
summary(fit)


#5.0 now getting into the actual comparison bit - looking at MassSA
#looked at percent mass lost, which is different, but standardizing for surface area fixes that
plot(massSA ~ taxon, data=pData)
points(massSA ~ taxon, data=pData)
mtext('Mass lost/surface area (mg/mm2)', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2) 

#percent lost by taxon by sealer presence
pData$seal <- factor(pData$seal)
noDef <- subset(pData, pData$taxon != 'Aragonite')
noDef <- subset(noDef, noDef$taxon != 'Calcite')
noDef <- subset(noDef, noDef$taxon != 'Centrostephanus')
droplevels(noDef$taxon)
boxplot(massSA ~ seal + droplevels(noDef$taxon), data=noDef, ylab = 'Mass/Surface Area', xlab = 'Taxon and Sealer')
stripchart(massSA ~ seal + droplevels(noDef$taxon), data = noDef, vertical=TRUE, add=TRUE, pch=1)


#6.0 Now run loop that does t test of paired seal/unsealed for each taxon
#both t test and power analysis in one loop to see if it is significant and if so, how much
ttestResult <- data.frame(taxon = c('Anadara','Pecten','Saccostrea') , p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd = NA, power = NA, shap=NA)
#note: 18 ana no, 11 ana yes, 18 pect no, 9 pect yes, 12 sacc no, 12 sacc yes
ttestResult$n <- 11
ttestResult[2, 'n'] <- 9
ttestResult[3, 'n'] <- 12

p = 0
for (T in ttestResult$taxon) {
  temp <- noDef[(noDef$taxon == T),]
  #temp$massSA = log(temp$massSA)
  ttest.temp <- t.test(massSA ~ seal, data=temp)
  shap.test <- shapiro.test(temp$massSA)
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=ttestResult[p, 'n'] ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult[p, 'mean.diff'] <- round(mean.diff, digits=3)
  #ttestResult[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  #ttestResult[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  ttestResult[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  ttestResult[p, 'a.p.value'] <- round(p.adjust(ttest.temp$p.value), digits=3)
  ttestResult[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult[p, 'sd'] <-  round(expSDBig[2], digits=4)
  ttestResult[p, 'shap'] <- round(shap.test$p.value, digits=3)
}
ttestResult$a.p.value <- round(p.adjust(ttestResult$p.value), digits=4)
print(ttestResult)

#7.0 Second loop that looks at the power of each of those t tests
#and how many samples would be needed to get power of > 0.8

ttestResult2 <- data.frame(taxon = c('Anadara','Pecten','Saccostrea') , p.value = NA, n = NA, power = 0.8)
p = 0
for (T in ttestResult2$taxon) {
  temp <- noDef[(noDef$taxon == T),]
  #temp$massSA = log(temp$massSA)
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$seal),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=NULL ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=0.8, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult2[p, 'p.value'] <- 0.05
  ttestResult2[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult2[p, 'n'] <-  power.temp$n
}
ttestResult2


#8.0 changing it into %/% standard
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

#now just the sealer, % standardized 
pData$seal <- factor(pData$seal)
noDef <- subset(pData, pData$taxon != 'Aragonite')
noDef <- subset(noDef, noDef$taxon != 'Calcite')
noDef <- subset(noDef, noDef$taxon != 'Centrostephanus')
droplevels(noDef$taxon)

#set colours using this line here
myCol <- ifelse(levels(factor(noDef$seal))=="Y" , rgb(0.1,0.1,0.7,0.5) , 
                ifelse(levels(factor(noDef$seal))=="N", rgb(0.8,0.1,0.3,0.6),
                       "grey90" ) )
pdf('test.pdf', width = 7, height=5)
boxplot(perMassSA ~ seal + droplevels(noDef$taxon), data=noDef, ylab = 'Standardized Mass/Surface Area', xlab = 'Taxon', col=myCol, xaxt = 'n', ylim = c(0.6,1.7))
stripchart(perMassSA ~ seal + droplevels(noDef$taxon), data = noDef, vertical=TRUE, add=TRUE, pch=1)
axis(1, at = c(1.5,3.5,5.5), labels=c(expression(italic('Anadara')), expression(italic('Pecten')), expression(italic('Saccostrea'))), las=1)
legend(x= 'topleft', legend = c('Unsealed','Sealed'), col= myCol, pch =15, cex=1.2)
dev.off()

#t-test stuff with % standardized
ttestResult <- data.frame(taxon = c('Anadara','Pecten','Saccostrea') , p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd = NA, power = NA, shap=NA)
ttestResult$n <- 11
ttestResult[2, 'n'] <- 9
ttestResult[3, 'n'] <- 12

p = 0
for (T in ttestResult$taxon) {
  temp <- noDef[(noDef$taxon == T),]
  #temp$perMassSA = log(temp$perMassSA)
  ttest.temp <- t.test(perMassSA ~ seal, data=temp)
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$seal),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$seal),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  shap.test <- shapiro.test(temp$perMassSA)
  p <- p + 1
  power.temp <- power.t.test(n=ttestResult[p, 'n'] ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult[p, 'mean.diff'] <- round(mean.diff, digits=3)
  #ttestResult[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  #ttestResult[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  ttestResult[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  ttestResult[p, 'a.p.value'] <- round(p.adjust(ttest.temp$p.value), digits=3)
  ttestResult[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult[p, 'sd'] <-  round(expSDBig[2], digits=4)
  ttestResult[p, 'shap'] <- round(shap.test$p.value, digits=3)
}
ttestResult$a.p.value <- round(p.adjust(ttestResult$p.value), digits=4)
ttestResult

#now how many with % standardized
ttestResult2 <- data.frame(taxon = c('Anadara','Pecten','Saccostrea') , p.value = NA, n = NA, power = 0.8)
p = 0
for (T in ttestResult2$taxon) {
  temp <- noDef[(noDef$taxon == T),]
  #temp$perMassSA = log(temp$perMassSA)
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$seal),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$seal),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=NULL ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=0.8, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult2[p, 'p.value'] <- 0.05
  ttestResult2[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult2[p, 'n'] <-  power.temp$n
}
ttestResult2

