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

#Script for creating graphics to visualize combined data from dissolution Master Dataset
#Portions of code adapted from 'Expt2 Post Processing Analysis'

#0.Open clean master dataset 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
taxaAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
TAXA2$taxon <- factor(TAXA2$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))
TAXA2 <- TAXA2[order(levels(TAXA2$taxon)),]
pData<-merge(pData,TAXA2,by='taxon')
pData$taxon <- factor(pData$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))


#1.1 Combining Expt 1 and Expt 2 without Expt 4
#first, you need to subset it
  
  Expt1_2 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T2.1'|pData$exptID == 'T2.2'|pData$exptID == 'T2.3')
  Expt1_2 <- subset(Expt123, Expt123$exptID == 'T1.1'|Expt123$exptID == 'T1.2'|Expt123$exptID == 'T1.3'|Expt123$exptID == 'T2.1'|Expt123$exptID == 'T2.2'|Expt123$exptID == 'T2.3')
  
#1.2 Setup and execution of split plots
  #pdf('./outFigs/1&2SplitWax.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  #1. firstly, some setup to subset
  library(ggplot2)
  factor(Expt1_2$waxYN)
  Expt1_2$waxYN <- as.factor(Expt1_2$waxYN)
  str(Expt1_2)
  
  ggplot(Expt1_2, aes(taxon, percentTotal/percentSATotal),) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic() + labs(x = 'Taxon') + 
    labs(y = '% Total Mass Lost/% Total Surface Area')
  
  dev.off()  
  
  
#2.1 Combining Expt 1 and Expt 4 without Expt 2, and then subsetting by experiment
  
  Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
  Expt1_3$exptNo <- 'Expt1'
  Expt1_3[(Expt1_3$exptID == 'T4.1'),'exptNo'] <-'Expt3'
  Expt1_3[(Expt1_3$exptID == 'T4.2'),'exptNo'] <-'Expt3'
  Expt1_3[(Expt1_3$exptID == 'T4.3'),'exptNo'] <-'Expt3'
  Expt1_3$massSA <- Expt1_3$cMass/Expt1_3$finalSA
 
  #setup %total 
  cTotal <- aggregate(Expt1_3$cMass, by = list(Expt1_3$exptID),FUN=sum)
  colnames(cTotal)<-c('exptID','cTotal')
  Expt1_3 <- merge(Expt1_3, cTotal, by = 'exptID')
  Expt1_3$percentTotal <- (Expt1_3$cMass)/(Expt1_3$cTotal) *100
  
  
  #now for surface area
  SATotal <- aggregate(Expt1_3$finalSA, by = list(Expt1_3$exptID),FUN=sum)
  colnames(SATotal)<-c('exptID','SATotal')
  Expt1_3 <- merge(Expt1_3, SATotal, by = 'exptID')
  Expt1_3$percentSATotal <- (Expt1_3$finalSA)/(Expt1_3$SATotal) * 100
  
  
#2.2 Setup and execution of comparison plots
  #only thing necessary for this is to split it into two more subsets by overall experiment 
  
  #pdf('./outFigs/1&3SplitExpt.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  
  #1. firstly, some setup to subset
  library(ggplot2)
  factor(Expt1_3$exptNo)
  Expt1_3$exptNo <- as.factor(Expt1_3$exptNo)
  str(Expt1_3)

  #then plot two side by side
  ggplot(Expt1_3, aes(taxon, massSA),) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic() + labs(x = 'Taxon') + 
    labs(y = 'Mass Lost/Surface Area (mg/mm\u00b2)')
  
  ggplot(Expt1_3, aes(taxon, percentTotal/percentSATotal),) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic() + labs(x = 'Taxon') + 
    labs(y = '% Total Mass Lost/% Total Surface Area')
 
#loop for running t tests and power analysis for all of the taxa: Mass/SA
  TAXA2 <- data.frame(taxon=taxa)
  summary1_3 <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd= NA, power = NA)
  summary1_3$n <- 15
    sub14 <- which(summary1_3$taxon == 'Marginopora'|summary1_3$taxon=='Turbo')
    summary1_3[sub14, 'n'] <-14
    sub13 <- which(summary1_3$taxon == 'Alaona')
    summary1_3[sub13, 'n'] <-13
    sub12 <- which(summary1_3$taxon == 'Liloa')
    summary1_3[sub12, 'n'] <-12
  
  p = 0
for (T in TAXA2$taxon) {
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  ttest.temp <- t.test(massSA~exptNo, data=temp)
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  power.temp <- power.t.test(n=15,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  p <- p + 1
  summary1_3[p, 'mean.diff'] <- round(mean.diff, digits=3)
  summary1_3[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  summary1_3[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  summary1_3[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  summary1_3[p, 'power'] <- round(power.temp$power, digits=2)
  summary1_3[p, 'sd'] <-  round(expSDBig[2], digits=2)
}

summary1_3$a.p.value <- round(p.adjust(summary1_3$p.value), digits=3)
print(summary1_3)


#next bit is trial that requires percentTotal/percent SA total calc  
Expt1_3$perMassSA <- Expt1_3$percentTotal/Expt1_3$percentSATotal
summary1_3Per <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , sample.size = NA, power = NA)
p = 0
for (T in TAXA2$taxon) {
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  ttest.temp <- t.test(perMassSA~exptNo, data=temp)
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=sd)
  #sample.temp <- ttest.temp$parameter + 1
  power.temp <- power.t.test(n=15,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDTemp[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  p <- p + 1
  summary1_3Per[p, 'mean.diff'] <- round(mean.diff, digits=3)
  summary1_3Per[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  summary1_3Per[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  summary1_3Per[p, 'p.value'] <- round(ttest.temp$p.value, digits=4)
  summary1_3Per[p, 'a.p.value'] <- round(p.adjust(ttest.temp$p.value), digits=4)
  summary1_3Per[p, 'power'] <- round(power.temp$power, digits=2)
}
summary1_3Per$a.p.value <- round(p.adjust(summary1_3Per$p.value), digits=4)
print(summary1_3Per)

#examining mass lost per surface area by taxon for expt1
  Expt1 <- subset(Expt123, Expt123$exptNo =='Expt1')
  Expt1$massSA <- Expt1$cMass/Expt1$finalSA
  plot(massSA~taxon, data=Expt1, ylab = 'Mass Lost/Surface Area (mg/mm\u00b3)', xlab = 'Taxon')
  points(massSA~taxon, data=Expt1)
  
  
#same for expt3
  Expt3 <- subset(Expt123, Expt123$exptNo =='Expt3')
  Expt3$massSA <- Expt3$cMass/Expt3$finalSA
  plot(massSA~taxon, data=Expt3, ylab = 'Mass Lost/Surface Area (mg/mm\u00b3)', xlab = 'Taxon')
  points(massSA~taxon, data=Expt3)
  modeltry2 <- kruskal.test(massSA~taxon, data = Expt3)
  p.adjust(modeltry2$p.value)
  