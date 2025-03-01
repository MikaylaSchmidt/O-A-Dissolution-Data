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
dShell <- read.csv('./Ch1Clean_MasterSet.csv')
#dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

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

#the natica fix + porosity  calc
#subNatNew <- which(pData$taxon == 'Notocochlis' & pData$waxYN == 'No Wax')
#pData[subNatNew,'finalSA']<- pData[subNatNew, 'finalSA'] * 1.5
#pData$poros <- pData$mass1/(pData$volume * 2.83)
#subCal <- which(pData$polymorph =='Calcite')
#pData[subCal, 'poros'] <- pData[subCal, 'mass1'] / (pData[subCal, 'volume'] * 2.711)
#pData$poros <- 1- pData$poros

#1.1 Combining Expt 1 and Expt 2 without Expt 4
#first, you need to subset it
  
  Expt1_2 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T2.1'|pData$exptID == 'T2.2'|pData$exptID == 'T2.3')
  #Expt1_2 <- subset(Expt123, Expt123$exptID == 'T1.1'|Expt123$exptID == 'T1.2'|Expt123$exptID == 'T1.3'|Expt123$exptID == 'T2.1'|Expt123$exptID == 'T2.2'|Expt123$exptID == 'T2.3')
  
#1.2 Setup and execution of split plots
  #pdf('./outFigs/1&2SplitWax.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  #1. firstly, some setup to subset
  library(ggplot2)
  factor(Expt1_2$waxYN)
  Expt1_2$waxYN <- as.factor(Expt1_2$waxYN)
  str(Expt1_2)
  
  #if you want to run this, need to run it LAST
  ggplot(Expt1_2, aes(taxon, percentTotal/percentSATotal),) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic() + labs(x = 'Taxon') + 
    labs(y = '% Total Mass Lost/% Total Surface Area')
  
  dev.off()  
  
  
#2.1 Combining Expt 1 and Expt 3 without Expt 2, and then subsetting by experiment
  
  Expt1_3 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T3.1'|pData$exptID == 'T3.2'|pData$exptID == 'T3.3')
  Expt1_3$exptNo <- 'Expt1'
  Expt1_3[(Expt1_3$exptID == 'T3.1'),'exptNo'] <-'Expt3'
  Expt1_3[(Expt1_3$exptID == 'T3.2'),'exptNo'] <-'Expt3'
  Expt1_3[(Expt1_3$exptID == 'T3.3'),'exptNo'] <-'Expt3'
  #Expt1_3$massSA <- Expt1_3$cMass/Expt1_3$finalSA
  
  #setup %total 
  #cTotal <- aggregate(Expt1_3$cMass, by = list(Expt1_3$exptID),FUN=sum)
  #colnames(cTotal)<-c('exptID','cTotal')
  #Expt1_3 <- merge(Expt1_3, cTotal, by = 'exptID')
  #Expt1_3$percentTotal <- (Expt1_3$cMass)/(Expt1_3$cTotal) *100
  
  
  #now for surface area
  #SATotal <- aggregate(Expt1_3$finalSA, by = list(Expt1_3$exptID),FUN=sum)
  #colnames(SATotal)<-c('exptID','SATotal')
  #Expt1_3 <- merge(Expt1_3, SATotal, by = 'exptID')
  #Expt1_3$percentSATotal <- (Expt1_3$finalSA)/(Expt1_3$SATotal) * 100
  
  
#2.2 Setup and execution of comparison plots
  #only thing necessary for this is to split it into two more subsets by overall experiment 
  
  #pdf('./outFigs/1&3SplitExpt.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  
  #1. firstly, some setup to subset
  library(ggplot2)
  library(ggthemes)
  factor(Expt1_3$exptNo)
  Expt1_3$exptNo <- as.factor(Expt1_3$exptNo)
  str(Expt1_3)

  #then plot two side by side
  ggplot(Expt1_3, aes(taxon, massSA),) +
    geom_boxplot(aes(fill = exptNo)) +
    theme_clean() + labs(x = 'Taxon') +
    labs(y = 'Mass Lost/Surface Area (mg/mm\u00b2)') +
    scale_fill_discrete(name = 'pH by time', labels =c('5.1 for 48 hours','7.1 for 168 hours'))
 
#print first graph for final
  pdf('massSAbyexpt.pdf', width = 13, height=5)
  
  #quick color vector
  myCol <- ifelse(levels(Expt1_3$exptNo)=="Expt1" , rgb(0.1,0.1,0.7,0.5) , 
                     ifelse(levels(Expt1_3$exptNo)=="Expt3", rgb(0.8,0.1,0.3,0.6),
                            "grey90" ) )
  
#the same plot but not using ggplot
  boxplot(massSA ~ exptNo * taxon, data = Expt1_3, xaxt = 'n', xlab = "", ylab = '', yaxt = 'n', col= myCol, xlim=c(1,22))
  stripchart(massSA ~ factor(exptNo)*taxon, data=Expt1_3, vertical=TRUE, add=TRUE, pch=1)
  mtext('Mass Lost/Surface Area (mg/mm\U00b2)', side=2, line=3)
  #Expt1_3$exptTax <- Expt1_3$exptNo.Expt1_3$taxon
  axis(2, las=1)
  mtext('Taxon', side = 1, line = 2)
  mtext(~italic('Ethalia'), side=1, line=0.5, at=1.5)
  mtext(~italic('Notocochlis'), side=1, line=0.5, at=3.5)
  mtext(~italic('Liloa'), side=1, line=0.5, at=5.5)
  mtext(~italic('Turbo'), side=1, line=0.5, at=7.5)
  mtext(~italic('Alaona'), side=1, line=0.5, at=9.5)
  mtext(~italic('Pinguitellina'), side=1, line=0.7, at=11.5)
  mtext(~italic('Fustiaria'), side=1, line=0.5, at=13.5)
  mtext(~italic('Halimeda'), side=1, line=0.5, at=15.5)
  mtext(~italic('Marginopora'), side=1, line=0.7, at=17.5)
  mtext('Aragonite', side=1, line=0.5, at=19.5)
  mtext('Calcite', side=1, line=0.5, at=21.5)
  legend("topleft", legend = c("pH 5.1 for 48 hours","pH 7.1 for 168 hours") , 
         col = c(rgb(0.1,0.1,0.7,0.5) , rgb(0.8,0.1,0.3,0.6)), bty = "o", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.05))

dev.off()           
#now do the same plot but for percentTotal/percentSATotal
pdf('perMassSAbyexpt.pdf', width = 13, height=5)
  Expt1_3$perMassSA <- Expt1_3$percentTotal/Expt1_3$percentSATotal 
  boxplot(perMassSA ~ (taxon * exptNo), data = Expt1_3, xlab = "", ylab = '', xaxt = 'n', yaxt="n", col= myCol, xlim=c(1,22))
  stripchart(perMassSA ~ taxon * exptNo, data=Expt1_3, vertical=TRUE, add=TRUE, pch=1)
  
  #points(massSA ~ exptNo * taxon, data= Expt1_3)
    mtext('Standardized Mass Lost/Surface Area', side=2, line=3)
    axis(2, las=1)
    #axis(1, las=1)
    mtext('Taxon', side = 1, line = 3)
    mtext(~italic('Ethalia'), side=1, line=0.5, at=1.5)
    mtext(~italic('Notocochlis'), side=1, line=0.5, at=3.6)
    mtext(~italic('Liloa'), side=1, line=0.5, at=5.6)
    mtext(~italic('Turbo'), side=1, line=0.5, at=7.5)
    mtext(~italic('Alaona'), side=1, line=0.5, at=9.6)
    mtext(~italic('Pinguitellina'), side=1, line=0.7, at=11.6)
    mtext(~italic('Fustiaria'), side=1, line=0.5, at=13.5)
    mtext(~italic('Halimeda'), side=1, line=0.5, at=15.5)
    mtext(~italic('Marginopora'), side=1, line=0.7, at=17.6)
    mtext('Aragonite', side=1, line=0.5, at=19.6)
    mtext('Calcite', side=1, line=0.5, at=21.5)
    abline(h=1)
    legend("topleft", legend = c("pH 5.1 for 48 hours","pH 7.1 for 168 hours") , 
         col = c(rgb(0.1,0.1,0.7,0.5) , rgb(0.8,0.1,0.3,0.6)), bty = "o", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.05))
 dev.off() 
  
  
  #now if you wanted to do this in ggplot you could do it like this
  ggplot(Expt1_3, aes(taxon, massSA)) +
    geom_boxplot() +
    theme_clean()
  
  ggplot(Expt1_3, aes(taxon, percentTotal/percentSATotal),) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_clean() + labs(x = 'Taxon') + 
    labs(y = '% Total Mass Lost/% Total Surface Area') +
    scale_fill_discrete(name = 'pH by time', labels =c('5.1 for 48 hours','7.1 for 168 hours'))
  
 

 #loop for running t tests and power analysis for all of the taxa: Mass/SA
  TAXA2 <- data.frame(taxon=taxa)
  summary1_3 <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd= NA, power = NA, shap= NA)
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
  shap.test <- shapiro.test(temp$massSA)
  ttest.temp <- t.test(massSA~exptNo, data=temp)
  expMeanTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=sd)
  expSDBig <- pmax.int(expSDTemp$x)
  p <- p + 1
  power.temp <- power.t.test(n=summary1_3[p, 'n'],delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power= NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  summary1_3[p, 'mean.diff'] <- round(mean.diff, digits=3)
  summary1_3[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  summary1_3[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  summary1_3[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
  summary1_3[p, 'power'] <- round(power.temp$power, digits=2)
  summary1_3[p, 'sd'] <-  round(expSDBig[2], digits=2)
  summary1_3[p, 'shap'] <- round(shap.test$p.value, digits=3)
}

summary1_3$a.p.value <- round(p.adjust(summary1_3$p.value), digits=3)
print(summary1_3)
summary1_3$taxon <- factor(summary1_3$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))

#next bit is trial that requires percentTotal/percent SA total calc  
Expt1_3$perMassSA <- Expt1_3$percentTotal/Expt1_3$percentSATotal
summary1_3Per <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = summary1_3$n, power = NA)
p = 0
for (T in TAXA2$taxon) {
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  ttest.temp <- t.test(perMassSA~exptNo, data=temp)
  expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
  expSDTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=sd)
  p <- p + 1
  power.temp <- power.t.test(n=summary1_3Per[p, 'n'],delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDTemp[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
  mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
  summary1_3Per[p, 'mean.diff'] <- round(mean.diff, digits=3)
  summary1_3Per[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
  summary1_3Per[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
  summary1_3Per[p, 'p.value'] <- round(ttest.temp$p.value, digits=4)
  summary1_3Per[p, 'power'] <- round(power.temp$power, digits=2)
}
summary1_3Per$a.p.value <- round(p.adjust(summary1_3Per$p.value), digits=4)
print(summary1_3Per)
summary1_3Per$taxon <- factor(summary1_3Per$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))

#how much difference is there in turbo dissolution
(ttest.temp$estimate[1] - ttest.temp$estimate[2])/ttest.temp$estimate[1]


#checking for normality in expt 1 mass/SA
  Expt1 <- subset(Expt1_3, Expt1_3$exptNo =='Expt1')
  summary1 <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA)
  p = 0
  for (T in TAXA2$taxon) {
    temp <- Expt1[(Expt1$taxon == T),]
    shap.test <- shapiro.test(temp$massSA)
    p <- p + 1
    summary1[p, 'p.value'] <- round(shap.test$p.value, digits=3)
  }
  summary1$a.p.value <- round(p.adjust(summary1$p.value), digits=4)
  
   
#same for expt3
  Expt3 <- subset(Expt1_3, Expt1_3$exptNo =='Expt3')
  summary3 <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA)
  p = 0
  for (T in TAXA2$taxon) {
    temp <- Expt3[(Expt1$taxon == T),]
    shap.test <- shapiro.test(temp$massSA)
    p <- p + 1
    summary3[p, 'p.value'] <- round(shap.test$p.value, digits=3)
  }
  summary3$a.p.value <- round(p.adjust(summary3$p.value), digits=4)
  
#then perMassSA  
  Expt1 <- subset(Expt1_3, Expt1_3$exptNo =='Expt1')
  summaryper1 <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA)
  p = 0
  for (T in TAXA2$taxon) {
    temp <- Expt1[(Expt1$taxon == T),]
    shap.test <- shapiro.test(temp$perMassSA)
    p <- p + 1
    summaryper1[p, 'p.value'] <- round(shap.test$p.value, digits=3)
  }
  summaryper1$a.p.value <- round(p.adjust(summaryper1$p.value), digits=4)
  
#examining metrics but with just calcite/marginopora?
onlyCalc <- subset(Expt1_3, Expt1_3$polymorph == 'Calcite')
droplevels(as.factor(onlyCalc$polymorph))
droplevels(as.factor(onlyCalc$taxon))
onlyCalc$taxon <- factor(onlyCalc$taxon, levels = c('Marginopora',  'Calcite'))
boxplot(massSA ~ exptNo * taxon, data = onlyCalc, col= myCol)

