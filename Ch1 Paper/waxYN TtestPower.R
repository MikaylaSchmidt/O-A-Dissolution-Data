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
##wax Y/N t test data looped
#to make the process shorter.
#intro code from fig gen 2
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt2.csv')
#dShell <- read.csv('./Ch1Clean_MasterSet.csv')
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
tAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), taxaAbrev = tAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
pData<-merge(pData,TAXA2,by='taxon')


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
waxData$waxYN <- factor(waxData$waxYN)
waxData <- subset(waxData, taxon != 'Calcite')

waxData$taxon <- c(droplevels(waxData$taxon))
waxData$taxon <- as.factor(waxData$taxon)
taxa <- sort(unique(waxData$taxon))
TAXA3 <- data.frame(taxon=taxa, shape = NA)
TAXA3$taxon <- factor(TAXA3$taxon, levels = c('Ethalia','Notocochlis',  'Alaona', 'Pinguitellina'))
TAXA3 <- TAXA3[order(levels(TAXA3$taxon)),]
waxData$taxon <- factor(waxData$taxon, levels = c('Ethalia','Notocochlis', 'Alaona', 'Pinguitellina'))
str(waxData)



#loop to run all of this easily
waxData$massSA <- waxData$cMass/waxData$finalSA
ttestResult <- data.frame(taxon = TAXA3$taxon, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd = NA, power = NA, power.p.value=NA)

ttestResult$n <- 14
sub15 <- which(ttestResult$taxon == 'Ethalia'|ttestResult$taxon=='Notocochlis')
ttestResult[sub15, 'n'] <-15
p = 0
for (T in TAXA3$taxon) {
  temp <- waxData[(waxData$taxon == T),]
  ttest.temp <- t.test(massSA ~ waxYN, data=temp)
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$waxYN),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$waxYN),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=ttestResult[p, 'n'] ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult[p, 'mean.diff'] <- round(mean.diff, digits=3)
  ttestResult[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  ttestResult[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  ttestResult[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  ttestResult[p, 'a.p.value'] <- round(p.adjust(ttest.temp$p.value), digits=3)
  ttestResult[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult[p, 'sd'] <-  round(expSDBig[2], digits=4)
}

ttestResult$a.p.value <- round(p.adjust(ttestResult$p.value), digits=4)

#power of adjusted p values
p=0
for (T in TAXA3$taxon) {
  temp <- waxData[(waxData$taxon == T),]
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$waxYN),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$waxYN),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.p.value.temp <- power.t.test(n=ttestResult[p, 'n'], delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=ttestResult[p, 'a.p.value'],power=NULL, type = c('two.sample'))
  ttestResult[p, 'power.p.value'] <- round(power.p.value.temp$power, digits=2)  
}

print(ttestResult)





#just looking at percentTotal/percentSA: does it negate differences?
waxData$perMassSA <- waxData$percentTotal/waxData$percentSATotal
ttestResult <- data.frame(taxon = TAXA3$taxon, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd = NA, power = NA, power.p.value=NA)

ttestResult$n <- 14
sub15 <- which(ttestResult$taxon == 'Ethalia'|ttestResult$taxon=='Notocochlis')
ttestResult[sub15, 'n'] <-15
p = 0
for (T in TAXA3$taxon) {
  temp <- waxData[(waxData$taxon == T),]
  ttest.temp <- t.test(perMassSA ~ waxYN, data=temp)
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$waxYN),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$waxYN),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=ttestResult[p, 'n'] ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  ttestResult[p, 'mean.diff'] <- round(mean.diff, digits=3)
  ttestResult[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  ttestResult[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  ttestResult[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  ttestResult[p, 'a.p.value'] <- round(p.adjust(ttest.temp$p.value), digits=3)
  ttestResult[p, 'power'] <-  round(power.temp$power, digits=2)
  ttestResult[p, 'sd'] <-  round(expSDBig[2], digits=2)
}

ttestResult$a.p.value <- round(p.adjust(ttestResult$p.value), digits=4)

#power of adjusted p values
p=0
for (T in TAXA3$taxon) {
  temp <- waxData[(waxData$taxon == T),]
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$waxYN),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$waxYN),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.p.value.temp <- power.t.test(n=ttestResult[p, 'n'], delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=ttestResult[p, 'a.p.value'],power=NULL, type = c('two.sample'))
  ttestResult[p, 'power.p.value'] <- round(power.p.value.temp$power, digits=2)  
}

print(ttestResult)


