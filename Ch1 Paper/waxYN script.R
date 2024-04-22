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

#intro code from fig gen 2
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
#dShell <- read.csv('./shellsData_Expt2.csv')
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

#run prev code if needed
#first separate out expt 2 then add relevant expt 1 stuff
waxData <- subset(Expt123, exptNo=='Expt2')
#additional text after this is for incorporating expt 1 into this, can skip
#newSub1 <- subset(Expt123, exptNo=='Expt1' & taxon=='Alaona')
#waxData <- rbind(waxData, newSub1)
#newSub2 <- subset(Expt123, exptNo=='Expt1' & taxon=='Pinguitellina')
#waxData <- rbind(waxData, newSub2)
#newSub3 <- subset(Expt123, exptNo=='Expt1' & taxon=='Ethalia')
#waxData <- rbind(waxData, newSub3)
#newSub4 <- subset(Expt123, exptNo=='Expt1' & taxon=='Notocochlis')
#waxData <- rbind(waxData, newSub4)

waxData$waxYN <- factor(waxData$waxYN)
waxData$taxon <- c(droplevels(waxData$taxon))
str(waxData)
levels(waxData$taxon) <- c(levels(waxData$taxon), 'Abranda2', 'Pingui2', 'Eth2', 'Nat2')
subAb <-  which(waxData$taxon == 'Alaona' & waxData$waxYN == 'Wax')
waxData[subAb,'taxon'] <- 'Abranda2'
subPin <-  which(waxData$taxon == 'Pinguitellina' & waxData$waxYN == 'Wax')
waxData[subPin,'taxon'] <- 'Pingui2'
subEth <-  which(waxData$taxon == 'Ethalia' & waxData$waxYN == 'Wax')
waxData[subEth,'taxon'] <- 'Eth2'
subNat <-  which(waxData$taxon == 'Notocochlis' & waxData$waxYN == 'Wax')
waxData[subNat,'taxon'] <- 'Nat2'

waxData$taxon <- factor(waxData$taxon, levels=c('Ethalia', 'Eth2', 'Notocochlis', 'Nat2', 'Alaona', 'Abranda2', 'Pinguitellina', 'Pingui2', 'Calcite'))
plot(pMass~taxon, data=waxData)
plot(percentTotal~taxon, data=waxData)
waxData$massSA <- waxData$cMass/waxData$finalSA
plot(massSA~taxon, data=waxData, ylab='Mass Lost/Surface Area (mg/mm\u00b2)' , xlab ='', xaxt ='n')
axis(side=1, at= c(waxData$taxon), labels = waxData$waxYN, tick = TRUE)
mtext('Ethalia', side=1, line=2.5, at=1.5)
mtext('Notocochlis', side=1, line=2.5, at=3.5)
mtext('Alaona', side=1, line=2.5, at=5.5)
mtext('Pinguitellina', side=1, line=2.5, at = 7.5)
mtext('Calcite', side=1, line=2.5, at =9)
mtext('Taxon', side =1, line =4.0, at =5)
points(massSA ~ taxon, data=waxData[!is.na(waxData$massSA),])
#abline(h=mean(waxData$massSA))
#kruskal.test(massSA~taxon, data=waxData)
#median(waxData$massSA)

subAb2 <- subset(waxData, waxData$taxon== 'Alaona'| waxData$taxon == 'Abranda2')
subAb2$taxon <- c(droplevels(subAb2$taxon))
expSDTemp <- aggregate(subAb2$massSA, by=list(subAb2$waxYN),FUN=sd)
#str(subAb2)
#t.test(pMass ~ taxon, data=subAb2)
#wilcox.test(percentTotal ~ taxon, data=subAb2, conf.int=TRUE)
abTTest <- t.test(massSA ~ taxon, data=subAb2)

subPin2 <- subset(waxData, waxData$taxon== 'Pinguitellina'| waxData$taxon == 'Pingui2')
subPin2$taxon <- c(droplevels(subPin2$taxon))
#t.test(percentTotal ~ taxon, data=subPin2)
#wilcox.test(percentTotal ~ taxon, data=subPin2, conf.int=TRUE)
pinTTest <- t.test(massSA ~ taxon, data=subPin2)

subEth2 <- subset(waxData, waxData$taxon== 'Ethalia'| waxData$taxon == 'Eth2')
subEth2$taxon <- c(droplevels(subEth2$taxon))
#t.test(percentTotal~taxon, data=subEth2)
#wilcox.test(percentTotal ~ taxon, data=subEth2, conf.int=TRUE)
ethTTest <- t.test(massSA ~ taxon, data=subEth2)

subNat2 <- subset(waxData, waxData$taxon== 'Notocochlis'| waxData$taxon == 'Nat2')
subNat2$taxon <- c(droplevels(subNat2$taxon))
#t.test(pMass~taxon, data=subNat2)
#wilcox.test(massSA ~ taxon, data=subNat2, conf.int=TRUE)
natTTest <- t.test(massSA ~ taxon, data=subNat2)


##this is for ttest results - can be with Expt 1 & 2 or just with Expt2 but still same
str(natTTest)
ttestResult <- data.frame(taxon = NA, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , sample.size = NA, power = NA)
ttestResult[1,] <- c('Alaona', abTTest$p.value, NA, abs(abTTest$estimate[1]-abTTest$estimate[2]), round(abTTest$conf.int[1], digits=3), round(abTTest$conf.int[2], digits=3), 14, NA)
ttestResult[2,] <- c('Pinguitellina', pinTTest$p.value, NA, NA, round(pinTTest$conf.int[1], digits=4), round(pinTTest$conf.int[2], digits=4), 14, NA)
ttestResult[3,] <- c('Ethalia', ethTTest$p.value, NA, NA, NA, NA, 15, NA)
ttestResult[4,] <- c('Notocochlis', natTTest$p.value, NA, NA, NA, NA, 15, NA)
ttestResult$a.p.value <- p.adjust(ttestResult$p.value)

print(ttestResult)
print(natTTest)

#power t.test for each taxon, starting with alaona
expMeanAb <- aggregate(subAb2$massSA, by=list(subAb2$taxon),FUN=mean, na.action=na.omit)
expMedianAb <- aggregate(subAb2$massSA, by=list(subAb2$taxon),FUN=median, na.action=na.omit)
expSDAb <- aggregate(subAb2$massSA, by=list(subAb2$taxon),FUN=sd)
power.t.test(n=14,delta=(expMeanAb[1,'x']-expMeanAb[2,'x']),sd=expSDAb[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
power.t.test(n=14,delta=(expMeanAb[1,'x']-expMeanAb[2,'x']),sd=expSDAb[2,'x'], sig.level=0.154,power=NULL, type = c('two.sample'))

#and pingui
expMeanPin <- aggregate(subPin2$massSA, by=list(subPin2$taxon),FUN=mean, na.action=na.omit)
expMedianPin <- aggregate(subPin2$massSA, by=list(subPin2$taxon),FUN=median, na.action=na.omit)
expSDPin <- aggregate(subPin2$massSA, by=list(subPin2$taxon),FUN=sd)
power.t.test(n=14,delta=(expMeanPin[1,'x']-expMeanPin[2,'x']),sd=expSDPin[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
power.t.test(n=14,delta=(expMeanPin[1,'x']-expMeanPin[2,'x']),sd=expSDPin[2,'x'], sig.level=1.0,power=NULL, type = c('two.sample'))


#then ethalia
expMeanEth <- aggregate(subEth2$massSA, by=list(subEth2$taxon),FUN=mean, na.action=na.omit)
expMedianEth <- aggregate(subEth2$massSA, by=list(subEth2$taxon),FUN=median, na.action=na.omit)
expSDEth <- aggregate(subEth2$massSA, by=list(subEth2$taxon),FUN=sd)
power.t.test(n=15,delta=(expMeanEth[1,'x']-expMeanEth[2,'x']),sd=expSDEth[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
power.t.test(n=15,delta=(expMeanEth[1,'x']-expMeanEth[2,'x']),sd=expSDEth[2,'x'], sig.level=1.00,power=NULL, type = c('two.sample'))

#finally Natica
expMeanNat <- aggregate(subNat2$massSA, by=list(subNat2$taxon),FUN=mean, na.action=na.omit)
expMedianNat <- aggregate(subNat2$massSA, by=list(subNat2$taxon),FUN=median, na.action=na.omit)
expSDNat <- aggregate(subNat2$massSA, by=list(subNat2$taxon),FUN=sd)
power.t.test(n=15,delta=(expMeanNat[1,'x']-expMeanNat[2,'x']),sd=expSDNat[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
expMeanNat[1,'x']/expMeanNat[2,'x']

