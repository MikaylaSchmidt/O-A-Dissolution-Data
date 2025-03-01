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
#Ch 2 Master Dataset Comparison
#0.0 Setup for everything
  setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
  expt1 <- read.csv('./Ch2Clean_Expt1.csv')
  expt2 <- read.csv('./Ch2Clean_Expt2.csv')
  expt2 <- expt2[,-19]
  expt3 <- read.csv('./Ch2Clean_Expt3.csv')
  masterSet <- rbind(expt1, expt2, expt3)
  masterSet <- masterSet[,-1]
  
  #porosity  calc
  masterSet$poros <- masterSet$mass1/(masterSet$volume * 2.83)
  subCal <- which(masterSet$polymorph =='Calcite')
  masterSet[subCal, 'poros'] <- masterSet[subCal, 'mass1'] / (masterSet[subCal, 'volume'] * 2.711)
  #masterSet$poros <- 1- masterSet$poros
  
  
  #adjustment of centrostephanus surface area
  #doesn't actually change the mg findings
  #centro SA hard to approximate + breakage is frequent
  #subCone <- which(masterSet$taxon == 'Centrostephanus')
  #masterSet[subCone,'finalSA'] <- masterSet[subCone,'finalSA'] * 0.75
  #masterSet$massSA <- masterSet$cMass/masterSet$finalSA
  
  
  #%/% calc
  cTotal <- aggregate(masterSet$cMass, by = list(masterSet$exptID),FUN=sum)
  colnames(cTotal)<-c('exptID','cTotal')
  masterSet <- merge(masterSet, cTotal, by = 'exptID')
  masterSet$percentTotal <- (masterSet$cMass)/(masterSet$cTotal) *100
  
  #now for surface area
  SATotal <- aggregate(masterSet$finalSA, by = list(masterSet$exptID),FUN=sum)
  colnames(SATotal)<-c('exptID','SATotal')
  masterSet <- merge(masterSet, SATotal, by = 'exptID')
  masterSet$percentSATotal <- (masterSet$finalSA)/(masterSet$SATotal) * 100
  
  #now standardized %/%
  masterSet$perMassSA <- masterSet$percentTotal/masterSet$percentSATotal

  #read in and merge hanna data
  hannaData <- read.csv('./Ch2 ExptData.csv')
  masterSet <- merge(masterSet, hannaData, by = 'exptID')
  masterSet$deltaPH <- masterSet$pH.99 - masterSet$pH.01
  masterSet$phHour <- masterSet$deltaPH/masterSet$duration

#quick save for later
write.csv(masterSet, 'Ch2Clean_MasterSet.csv', row.names=TRUE)
#dev.off()

#quick look at spec type averages
subAna <- subset(masterSet, masterSet$taxon == 'Anadara')
hist(subAna$finalSA)
mean(subAna$finalSA)
mean(subAna$mass1)
subSac <- subset(masterSet, masterSet$taxon == 'Saccostrea')
hist(subSac$finalSA)
mean(subSac$finalSA)
mean(subSac$mass1)
subPect <- subset(masterSet, masterSet$taxon == 'Pecten')
hist(subPect$finalSA)
mean(subPect$finalSA)
mean(subPect$mass1)


#fixing issues with factoring
TAXA <- unique(masterSet$taxon)
masterSet$taxon <- as.factor(masterSet$taxon)
taxa <- sort(unique(masterSet$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Anadara')] <- 1
taxaAbrev <- substring(taxa,0,5)
TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
masterSet<-merge(masterSet,TAXA2,by='taxon')

  #set up expt no variable
  masterSet$exptNo <- 'Expt1'
  masterSet[(masterSet$exptID == '2.1'),'exptNo'] <-'Expt2'
  masterSet[(masterSet$exptID == '2.2'),'exptNo'] <-'Expt2'
  masterSet[(masterSet$exptID == '2.3'),'exptNo'] <-'Expt2'
  masterSet[(masterSet$exptID == '3.1'),'exptNo'] <-'Expt3'
  masterSet[(masterSet$exptID == '3.2'),'exptNo'] <-'Expt3'
  masterSet[(masterSet$exptID == '3.3'),'exptNo'] <-'Expt3'
  Expt123 <- masterSet
  
  
#1.0 Set up experimental variable by getting rid of expt 2  
 
  
  plot(massSA~taxon, data=masterSet)
  points(massSA~taxon, data=masterSet)
  factor(masterSet$exptNo)
  
  masterSet2 <- subset(masterSet, exptNo != 'Expt2')
  boxplot(massSA~masterSet2$exptNo+taxon, data=masterSet2)
  stripchart(massSA ~ masterSet2$exptNo+droplevels(taxon), data = masterSet2, vertical=TRUE, add=TRUE, pch=1)
  str(masterSet)
  
  plot(perMassSA~taxon, data=masterSet2)
  points(perMassSA~taxon, data=masterSet2)
  boxplot(perMassSA~factor(masterSet2$exptNo)+taxon, data=masterSet2)
  stripchart(perMassSA ~ factor(masterSet2$exptNo)+droplevels(taxon), data = masterSet2, vertical=TRUE, add=TRUE, pch=1)
  
  plot(massSA~taxon, data=masterSet2)
  points(massSA~taxon, data=masterSet2)
  boxplot(massSA~factor(masterSet2$exptNo)+taxon, data=masterSet2)
  stripchart(perMassSA ~ factor(masterSet2$exptNo)+droplevels(taxon), data = masterSet2, vertical=TRUE, add=TRUE, pch=1)
  
  
#2.0 looking to see if there's difference between expt 1 and 3
  Expt1_3 <- subset(masterSet, exptNo != 'Expt2')
  Expt1_3 <- subset(Expt1_3, taxon != 'Coral')
  Expt1_3 <- subset(Expt1_3, taxon != 'Marginopora')
  boxplot(massSA~factor(Expt1_3$exptNo)+droplevels(taxon), data=Expt1_3)
  stripchart(massSA ~ factor(Expt1_3$exptNo)+droplevels(taxon), data = Expt1_3, vertical=TRUE, add=TRUE, pch=1)
  
  #loop for running t tests and power analysis for all of the taxa: Mass/SA
  taxa2 <- sort(unique(droplevels(Expt1_3$taxon)))
  TAXA2 <- data.frame(taxon=taxa2)
  summary1_3 <- data.frame(taxon = taxa2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = NA, sd= NA, power = NA, shap= NA)
  summary1_3$n <- 17
  sub15 <- which(summary1_3$taxon == 'Saccostrea')
  summary1_3[sub15, 'n'] <-15
  sub7 <- which(summary1_3$taxon == 'Centrostephanus')
  summary1_3[sub7, 'n'] <-7
  
  p = 0
  for (T in TAXA2$taxon) {
    temp <- Expt1_3[(Expt1_3$taxon == T),]
    shap.test <- shapiro.test(temp$massSA)
    ttest.temp <- wilcox.test(massSA~exptNo, data=temp)
    expMeanTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
    expSDTemp <- aggregate(temp$massSA, by=list(temp$exptNo),FUN=sd)
    expSDBig <- pmax.int(expSDTemp$x)
    p <- p + 1
    #power.temp <- power.t.test(n=summary1_3[p, 'n'],delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power= NULL, type = c('two.sample'))
    mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
    summary1_3[p, 'mean.diff'] <- round(mean.diff, digits=3)
    #summary1_3[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
    #summary1_3[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
    summary1_3[p, 'p.value'] <- round(ttest.temp$p.value, digits=3)
    #summary1_3[p, 'power'] <- round(power.temp$power, digits=2)
    summary1_3[p, 'sd'] <-  round(expSDBig[2], digits=2)
    summary1_3[p, 'shap'] <- round(shap.test$p.value, digits=3)
  }
  
  summary1_3$a.p.value <- round(p.adjust(summary1_3$p.value), digits=3)
  print(summary1_3)
  
  
#3.0 Now onto %/% standard
  
  #set colours using this line here
  myCol <- ifelse(levels(factor(Expt1_3$exptNo))=="Expt1" , rgb(0.1,0.1,0.7,0.5) , 
                  ifelse(levels(factor(Expt1_3$exptNo))=="Expt3", rgb(0.8,0.1,0.3,0.6),
                         "grey90" ) )
  
  #boxplot of the things you actually want to look at
  boxplot(perMassSA~exptNo * droplevels(taxon), data= Expt1_3, ylab = 'Mass/Surface Area (Standardized)', xlab = 'Taxon', col= myCol, xaxt ='n')
  axis(1, at = c(1.5,3.5,5.5,7.5,9.5,11.5), labels=c(expression(italic('Anadara')), 'Aragonite', 'Calcite', expression(italic('Centrostephanus')), expression(italic('Pecten')), expression(italic('Saccostrea'))), las=1)
  stripchart(perMassSA ~ factor(Expt1_3$exptNo)+droplevels(taxon), data = Expt1_3, vertical=TRUE, add=TRUE, pch=1)
  legend(x= 'topleft', legend = c('Expt1','Expt3'), col= myCol, pch = 19)
  abline(h=1)
  
  #next bit is trial that requires percentTotal/percent SA total calc  
  Expt1_3$perMassSA <- Expt1_3$percentTotal/Expt1_3$percentSATotal
  summary1_3Per <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = 17, SD = NA, power = NA, shap = NA)
  sub15 <- which(summary1_3Per$taxon == 'Saccostrea')
  summary1_3Per[sub15, 'n'] <-15
  sub7 <- which(summary1_3Per$taxon == 'Centrostephanus')
  summary1_3Per[sub7, 'n'] <-7
  
  p = 0
  for (T in TAXA2$taxon) {
    temp <- Expt1_3[(Expt1_3$taxon == T),]
    ttest.temp <- wilcox.test(perMassSA~exptNo, data=temp)
    shap.test <- shapiro.test(temp$perMassSA)
    expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
    expSDTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=sd)
    p <- p + 1
    power.temp <- power.t.test(n=summary1_3Per[p, 'n'],delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDTemp[2,'x'], sig.level=0.05,power=NULL, type = c('two.sample'))
    mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
    summary1_3Per[p, 'mean.diff'] <- round(mean.diff, digits=3)
    #summary1_3Per[p, 'conf.int1'] <- round(ttest.temp$conf.int[1], digits=3)
    #summary1_3Per[p, 'conf.int2'] <- round(ttest.temp$conf.int[2], digits=3)
    #summary1_3Per[p, 'SD'] <- round(ttest.temp$stderr, digits=4)
    summary1_3Per[p, 'p.value'] <- round(ttest.temp$p.value, digits=4)
    #summary1_3Per[p, 'power'] <- round(power.temp$power, digits=2)
    summary1_3Per[p, 'shap'] <- round(shap.test$p.value, digits=3)
  }
  summary1_3Per$a.p.value <- round(p.adjust(summary1_3Per$p.value), digits=4)
  print(summary1_3Per)
  
  #summary1_3Per$taxon <- factor(summary1_3Per$taxon, levels = c('Ethalia','Notocochlis', 'Liloa', 'Turbo', 'Alaona', 'Pinguitellina', 'Fustiaria', 'Halimeda', 'Marginopora', 'Aragonite', 'Calcite'))

#checking for normality for this info can be found in the Master Dataset Comparison file for CH1
  #if you need to do that
  
#second summary1_3 where you look at how many specs would be needed to tell difference 
  summary1_3Per2 <- data.frame(taxon = TAXA2 , p.value = NA, n = NA, power = 0.8)
  p = 0
  for (T in summary1_3Per2$taxon) {
    temp <- Expt1_3[(Expt1_3$taxon == T),]
    expMeanTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=mean, na.action=na.omit)
    expSDTemp <- aggregate(temp$perMassSA, by=list(temp$exptNo),FUN=sd)
    expSDBig <- pmax.int(expSDTemp$x)
    p <- p + 1
    power.temp <- power.t.test(n=NULL ,delta=(expMeanTemp[1,'x']-expMeanTemp[2,'x']),sd=expSDBig[2], sig.level=0.05,power=0.8, type = c('two.sample'))
    mean.diff <- abs(expMeanTemp[1,'x']-expMeanTemp[2,'x'])
    summary1_3Per2[p, 'p.value'] <- 0.05
    summary1_3Per2[p, 'power'] <-  round(power.temp$power, digits=2)
    summary1_3Per2[p, 'n'] <-  power.temp$n
  }
  summary1_3Per2
  
#4.0  Matrix of which ones are actually different compared to each other, not separated by taxon
#First, code to generate the visual of what you're looking at: which of these are statistically different?
  #set colours using this line here
 
  masterSet$taxon <- factor(masterSet$taxon, levels = c('Anadara','Pecten', 'Saccostrea', 'Centrostephanus', 'Marginopora', 'Coral', 'Aragonite', 'Calcite'))
  
  
  pdf('ch2taxon.pdf', width = 12, height=7)
  plot(perMassSA~taxon, data=masterSet, xlab = 'Taxon', ylab = 'Standardized Mass Lost/Surface Area', xaxt = 'n', xlim=c(0.5,8.5))
  #axis(1,at=1:length(taxa), labels = taxa)
  axis(1, at = c(1,2,3,4,5,6,7,8),labels =c('','','','','','','',''))
       #labels=c(expression(italic('Anadara')), 'Aragonite', 'Calcite', expression(italic('Centrostephanus')),expression(italic('Seriatopora')), expression(italic('Marginopora')), expression(italic('Pecten')), expression(italic('Saccostrea'))), las=1)
  mtext(~italic('Anadara'), side=1, line=0.75, at=1)
  mtext(~italic('Pecten'), side=1, line=0.75, at=2)
  mtext(~italic('Saccostrea'), side=1, line=0.75, at=3)
  mtext(~italic('Centrostephanus'), side=1, line=0.85, at=4)
  mtext(~italic('Marginopora'), side=1, line=0.85, at=5)
  mtext(~italic('Seriatopora'), side=1, line=0.85, at=6)
  mtext('Aragonite', side=1, line=0.7, at=7)
  mtext('Calcite', side=1, line=0.7, at=8)
  points(perMassSA~taxon, data=masterSet)
  abline(h=1)
dev.off()
  
  plot(pMass~taxon, data=masterSet, xlab = 'Taxon', ylab = 'Percent Mass Lost')
  points(pMass~taxon, data=masterSet)
  
  #matrix for running t tests between all of the different specs
master.frame <- data.frame(taxon= masterSet$taxon, perMassSA= masterSet$perMassSA)
expMean <- aggregate(masterSet$perMassSA, by=list(masterSet$taxon),FUN=mean, na.action=na.omit)

  tax <- vector()
  res <- vector()
  c <- 0
  
  for (a in 1:nrow(expMean)) {
    for (b in a:nrow(expMean)) {
      
      c <- c + 1
      tax[c] <- paste(expMean[a,'Group.1'],expMean[b,'Group.1'])
      
      if (a == b) {
        res[c] <- -99
      } else {
        temp <- master.frame[(master.frame$taxon == expMean[a, 'Group.1']),]
        temp2 <- master.frame[(master.frame$taxon == expMean[b, 'Group.1']),]
        temp3 <- rbind(temp,temp2)
        #ttest.temp <- t.test(perMassSA ~ taxon, data = temp3)
        ttest.temp <- wilcox.test(perMassSA ~ taxon, data = temp3)
        res[c] <- ttest.temp$p.value
      }
    }}
  
  dataTrial <- data.frame(comparison=tax, p.value=round(res, digits=6), adj.p.value =NA)
  dataTrial <- dataTrial[(dataTrial$p.value != -99),]
  dataTrial$adj.p.value <- p.adjust(dataTrial$p.value)
  dataTrial
  dataTrialSig <- dataTrial[(dataTrial$adj.p.value < 0.05),]
  
##compare all to 1
  
  
#5.0 Looking at power of each ttest?
  
  #all of means and things needed for matrix of power t tests
  expMean <- aggregate(masterSet$perMassSA, by=list(masterSet$taxon),FUN=mean, na.action=na.omit)
  expSD <- aggregate(masterSet$perMassSA, by=list(masterSet$taxon),FUN=sd)
  
  #longer script to run power analysis of t test for number of specimens needed to discern massSA differences if present
  tax <- vector()
  res <- vector()
  c <- 0
  
  for (a in 1:nrow(expMean)) {
    for (b in a:nrow(expMean)) {
      
      c <- c + 1
      tax[c] <- paste(expMean[a,'Group.1'],expMean[b,'Group.1'])
      
      if (a == b) {
        res[c] <- -99
      } else {
        
        DELTA <- expMean[a,'x'] - expMean[(b),'x']
        SD <- expSD[a,'x']
        if (expSD[(b),'x'] > SD)
          SD <- expSD[(b),'x']
        
        res[c] <- power.t.test(n=10,delta=DELTA,sd=SD, sig.level=0.05,power=NULL)$power
      }
    }}
  
  dataTrial <- data.frame(comparison=tax, power=res)
  dataTrial <- dataTrial[(dataTrial$power != -99),]
  dataTrial
  

#6.0 Comparing stuff - what if you look at just calcite specs? Does that improve things???????
#first, calcite subset
  onlyCalc <- subset(masterSet, masterSet$polymorph == 'Calcite')
  droplevels(as.factor(onlyCalc$polymorph))
  droplevels(as.factor(onlyCalc$taxon))
  onlyCalc$taxon <- factor(onlyCalc$taxon, levels = c('Pecten',  'Saccostrea', 'Centrostephanus', 'Marginopora', 'Calcite'))
  
#before modeling, a visual of the %/% mass stuff
  plot(log2(percentTotal) ~ log2(percentSATotal), data=onlyCalc, col=tColor, pch=tPoint)
  line1 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=onlyCalc)
  abline(line1)
  summary(line1)
  
  
  
#load in library stuff  
library(lmerTest)
library(MuMIn)
library(cAIC4)

 
  
 #6.1 plot of massSA
  boxplot(massSA ~ taxon, data = onlyCalc)
  points(massSA~taxon, data=onlyCalc)
  
  #is mass/SA plot better fit by only calcite? 
  plot(cMass ~ finalSA, data=onlyCalc)
  only1 <- lm(cMass ~ finalSA, data=onlyCalc)
  summary(only1)
  plot(log2(cMass) ~ log2(finalSA), data=onlyCalc)
  only1.1 <- lm(log2(cMass) ~ log2(finalSA), data=onlyCalc)
  summary(only1.1)
  only1.2 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|exptID) + (1|taxon) , data = masterSet)
  summary(only1.2)
  r.squaredGLMM(only1.2)
  only1.3 <- lm(log2(cMass) ~ log2(deviation), data=onlyCalc)
  summary(only1.3)
  
  #plot of perMassSA
  boxplot(perMassSA ~taxon, data = onlyCalc, xlab='Taxon', ylab ="Mass/Surface Area (Standardized)")
  points(perMassSA~taxon, data=onlyCalc)
  
  #perMassSA plot fit with just calcite 
  #this is bad - if you look at plot of linear model, data is skewed
  #so you should log(2) it. see only3 plot, is much better than only2
  
  plot(percentTotal ~ percentSATotal, data=onlyCalc, ylab = 'Standardized Mass', xlab='Standardized Surface Area')
  only2 <- lm(percentTotal ~ percentSATotal, data=onlyCalc)
  summary(only2)
  plot(only2)
  
  only3 <- lm(log2(percentTotal) ~ log2(percentSATotal), data=onlyCalc)
  summary(only3)
  plot(only3)
  
  #the addition of magnesium to the equation. if you remove saccostrea it is very important
  only4 <- lm(log2(percentTotal) ~ log2(percentSATotal) + log2(mag), data=onlyCalc)
  summary(only4)

  subNoSacc <- subset(onlyCalc, taxon !='Saccostrea')
  only5 <- lm(log2(percentTotal) ~ log2(percentSATotal) + log2(mag), data=subNoSacc)
  summary(only5)
  
#6.2 create data frame with only the relevant variables from onlyCalc 
#first modelling using the %/% method
  relDataCalc <- onlyCalc[,c('cMass', 'pMass', 'mass1', "mass2", 'percentTotal', 'finalSA', 'percentSATotal', 'volume', 'cSize', 'poros', 'deviation', 'mag', 'taxon', 'exptID', 'pH.01', 'pH.99', 'meanTemp', 'deltaPH', 'duration', 'phHour')]
  
  cal1 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|exptID) + (1|taxon) , data = relDataCalc)
  summary(cal1)
  rand(cal1)
  r.squaredGLMM(cal1)
  
  cal2 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|exptID), data = relDataCalc)
  summary(cal2)
  rand(cal2)
  r.squaredGLMM(cal2)

  cal2.1 <- lmer(log2(cMass) ~  log2(finalSA) + log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|taxon), data = relDataCalc)
  summary(cal2.1)
  rand(cal2.1)
  r.squaredGLMM(cal2.1)
  
  cal2.2 <- lm(log2(cMass) ~  log2(finalSA) + log2(volume) + poros + log2(deviation) + log2(mag), data = relDataCalc)
  dredge(cal2.2)
  summary(cal2.2)
  rand(cal2.2)
  r.squaredGLMM(cal2.2)
  
  
#quick correlation matrix
  allNum1 <-relDataCalc[,c('percentSATotal', 'volume', 'cSize', 'poros', 'deviation', 'mag', 'pH.01', 'meanTemp', 'deltaPH')]
  cor1 <- cor(allNum1)
  print(cor1)
  xtable(cor1)
  cor1 <- cor1[(cor1 > 0.70)]
  print(cor1)
  allNum2 <-relDataCalc[,c('pMass', 'finalSA', 'volume', 'cSize', 'poros', 'deviation', 'mag', 'pH.01', 'meanTemp', 'deltaPH')]
  cor2 <- cor(allNum2)
  print(cor2)
  cor2 <- cor2[(cor2 > 0.70)]
  print(cor2)
  
#6.3 cal3 is the kitchen sink one. nothing is significant except for meanTemp. exptID not sig if you have meanTemp 
  ##for some reason meanTemp causes the intercept to shit itself
  cal3 <- lmer(log2(percentTotal) ~  log2(percentSATotal)  + log2(volume) + log2(poros) + log2(deviation) + log2(mag)+ log2(meanTemp) + log2(deltaPH) + log2(pH.01) + (1|exptID) + (1|taxon), data = relDataCalc)
  summary(cal3)
  rand(cal3)
  r.squaredGLMM(cal3)
  
  #now for beta coefficients for full model 
  relDataCalcScale <- relDataCalc
  relDataCalcScale[,c(5:11,12,15:20)] <- scale(log2(relDataCalc[,c(5:11,12,15:20)]))
  
  cal3beta <- lmer(percentTotal ~  percentSATotal +  volume + poros + deviation + mag + meanTemp + deltaPH + pH.01 + (1|taxon) + (1|exptID), data = relDataCalcScale)
  summary(cal3beta)
  round((summary(cal3beta)$coefficients)^2, digits=5)
  
#6.4 exptID only has 1% difference on r2 and soaks up the effect of meanTemp
  #this also shows how variable for taxon absorbs stuff from volume, poros, and deviation
  #below is with meantemp
  cal3.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + log2(poros) + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|exptID), data = relDataCalc)
  summary(cal3.1)
  rand(cal3.1)
  r.squaredGLMM(cal3.1)

  #just the lm. make comparable by computing beta coefficients
  cal3.1Lin <- lm(log2(percentTotal) ~  log2(percentSATotal) + log2(volume) + log2(poros) + log2(deviation) + log2(mag) + log2(meanTemp) + log2(deltaPH) + log2(pH.01), data = relDataCalc)
  summary(cal3.1Lin)
  options(na.action = "na.fail")
  dredge(cal3.1Lin)
  
  cal3.1beta <- lm(percentTotal ~  percentSATotal + volume + poros + deviation + mag + meanTemp + deltaPH + pH.01, data = relDataCalcScale)
  summary(cal3.1beta)
  round((summary(cal3.1beta)$coefficients)^2, digits=3)
  #best model
  calBest <- lm(log2(percentTotal) ~  log2(percentSATotal) + log2(poros) + log2(deviation), data = relDataCalc)
  summary(calBest)
  plot(calBest)
  calBestbeta <- lm(percentTotal ~  percentSATotal + poros + deviation, data = relDataCalcScale)
  summary(calBestbeta)

  #and then without expt variables
  cal3.2 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + log2(mag) + (1|exptID), data = relDataCalc)
  summary(cal3.2)
  rand(cal3.2)
  r.squaredGLMM(cal3.2)
  
  
#6.5 kitchen sink model minus exptID, to look at just taxon
  #meanTemp is postive as expected
  cal4 <- lmer(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag) + (1|taxon), data = relDataCalc)
  summary(cal4)
  rand(cal4)
  r.squaredGLMM(cal4)
  
  #then without measured variables other than percentSATotal
  cal4.1 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + (1|taxon), data = relDataCalc)
  summary(cal4.1)
  rand(cal4.1)
  r.squaredGLMM(cal4.1)
  
  #just taxon.
  cal4.2 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + (1|taxon), data = relDataCalc)
  summary(cal4.2)
  rand(cal4.2)
  r.squaredGLMM(cal4.2)

  #just exptID.
  cal5 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + (1|exptID), data = relDataCalc)
  summary(cal5)
  r.squaredGLMM(cal5)

#6.6 final model: quick lm with all things other than taxon. effects from: percentSATotal, poros, deviation, meanTemp, and mag and with an r2 of 85%. means taxon soaks up poros, deviation, and magnesium plus an extra 2%
  cal5.1 <- lm(log2(percentTotal) ~  log2(percentSATotal) +  log2(cSize) + log2(volume) + poros + log2(deviation) + log2(meanTemp) + log2(deltaPH) + log2(duration) + log2(pH.01) + log2(mag), data = relDataCalc)
  summary(cal5.1)

#if you pare things down to only the significant variables, mag is significant 
  cal5.2 <- lm(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation) + log2(mag), data = relDataCalc)
  summary(cal5.2)
  
#this is where things get confusing
  cal5.3 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation) + log2(mag) + (1|exptID) + (1|taxon), data = relDataCalc)
  summary(cal5.3)
  rand(cal5.3)
  r.squaredGLMM(cal5.3)
  
  
  cal5.4 <- lmer(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation) + log2(mag) + (1|taxon), data = relDataCalc)
  summary(cal5.4)
  rand(cal5.4)
  r.squaredGLMM(cal5.4)
  
#6.7 quick lm 
  cal5.5 <- lm(log2(percentTotal) ~  log2(percentSATotal) + poros + log2(deviation), data = relDataCalc)
  summary(cal5.5)
  
  cal5.6 <- lm(log2(percentTotal) ~  log2(percentSATotal), data = relDataCalc)
  summary(cal5.6)
  
  
  
#7.0 test using the john model of pMass lost

  cal6 <- lmer(log2(mass2/mass1) ~  log2(finalSA) + log2(cSize) + poros + log2 (deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(pH.99) + log2(mag) + (1|exptID) + (1|taxon), data = relDataCalc)
  summary(cal6)
  r.squaredGLMM(cal6)
  
  cal7 <- lmer(log2(pMass) ~  log2(finalSA) + log2(cSize) + poros + log2 (deviation) + log2(volume) + log2(meanTemp) + log2(deltaPH) + log2(pH.01) + log2(pH.99) + log2(mag) + (1|exptID) + (1|taxon), data = relDataCalc)
  summary(cal7)
  r.squaredGLMM(cal7)
  
  #comparing the two reverse ones to see what is normally distributed 
  plot(log2(mass2/mass1) ~ log2(pMass), data=onlyCalc)
  plot(mass2/mass1 ~pMass, data=onlyCalc)
  hist(1-onlyCalc$pMass)
  hist(onlyCalc$mass2/onlyCalc$mass1)
  hist(log2(onlyCalc$pMass))
  hist(log2(onlyCalc$mass2/onlyCalc$mass1))
  
  #what other ones are normally distributed? only poros but might as well log it
  hist(onlyCalc$pMass)
  hist(onlyCalc$finalSA)
  hist(onlyCalc$volume)
  hist(onlyCalc$deviation)
  hist(onlyCalc$poros)
  
  testModel<- lm(log2(pMass)~log2(finalSA), data=onlyCalc)
  summary(testModel)
  
  testModel2<- lm(log2(pMass)~log2(finalSA) + log2(volume), data=onlyCalc)
  summary(testModel2)
  

  testModel3<- lmer(log2(pMass)~log2(finalSA) + log2(cSize) + (1|taxon), data=onlyCalc)
  summary(testModel3)
  ranef(testModel3)
  rand(testModel3)
  r.squaredGLMM(testModel3)

#7.1 full kitchen sink model here. three things important: finalSA, volume, and porosity. also, exptID and taxon have effects to the tune of 15% of variance
  testModel4<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(poros) + log2(deviation) + log2(mag) + log2(meanTemp)+ log2(pH.01) + log2(deltaPH)  + (1|exptID) + (1|taxon) , data=onlyCalc)
  summary(testModel4)
  rand(testModel4)
  r.squaredGLMM(testModel4)
  
  #setup for beta coefficients
  relDataCalcScale <- relDataCalc
  relDataCalcScale[,c(2,5:9,11,12,15:20)] <- scale(log2(relDataCalc[,c(2,5:9,11,12,15:20)]))
  relDataCalcScale$poros <- scale(relDataCalcScale$poros)
  testModel4beta <- lmer(pMass ~  finalSA + volume + poros + deviation + mag + meanTemp + deltaPH + pH.01 + (1|taxon) + (1|exptID), data = relDataCalcScale)
  summary(testModel4beta)
  round((summary(testModel4beta)$coefficients)^2, digits=3)
  
  #same model but minus random effects, in lm form
  testModel4.1<- lm(log2(pMass)~log2(finalSA) + log2(volume) +  poros + log2(deviation) + log2(mag) + log2(meanTemp)+ log2(pH.01) + log2(deltaPH), data=onlyCalc)
  summary(testModel4.1)
  testModel4.1beta <- lm(pMass ~  finalSA + volume + poros + deviation + mag + meanTemp + deltaPH + pH.01, data = relDataCalcScale)
  summary(testModel4.1beta)
  round((summary(testModel4.1beta)$coefficients)^2, digits=3)
  
  
#7.2 below shows that expt ID is significant but only changes things by like 4% in the r2
  testModel5<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(mag)+ log2(deltaPH) + log2(meanTemp) + log2(pH.01) + (1|exptID), data=onlyCalc)
  summary(testModel5)
  rand(testModel5)
  r.squaredGLMM(testModel5)
  
  #exptID yes, expt variables no
  testModel5.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(mag) + (1|exptID), data=onlyCalc)
  summary(testModel5.1)
  rand(testModel5.1)
  r.squaredGLMM(testModel5.1)
  
#7.3 exptID also soaks up the meanTemp + pH.01 difference issues seen in non exptID when you don't include it
  #taxon is very important though. without taxon effect this only accounts for 79% of variance vs 88% with taxon
  #main predictors otherwise are surface area, volume, and porosity
  #also indicates that volume and cSize are highly correlated, which makes sense
  testModel6<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + log2(cSize) + log2(deviation) + poros + log2(mag)+ log2(deltaPH) + log2(meanTemp) + log2(pH.01) + (1|taxon), data=onlyCalc)
  summary(testModel6)
  rand(testModel6)
  r.squaredGLMM(testModel6)
  
  testModel6.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + poros + log2(mag)+ log2(deltaPH) + log2(meanTemp) + log2(pH.01) + (1|taxon), data=onlyCalc)
  summary(testModel6.1)
  rand(testModel6.1)
  r.squaredGLMM(testModel6.1)

#7.4 finally, a model without taxon. this is probably the ideal model. mag not significant.
#why does this account for so much rsquared compared to the r.squared glmm from testModel6???
  
  testModelIdeal<- lm(log2(pMass)~log2(finalSA) + log2(volume) + log2(poros), data=onlyCalc)
  summary(testModelIdeal)
  
  testModel8.1<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + poros + (1|taxon), data=onlyCalc)
  summary(testModel8.1)
  rand(testModel8.1)
  r.squaredGLMM(testModel8.1)
  
  testModel8.2<- lmer(log2(pMass)~log2(finalSA) + log2(volume) + poros + (1|taxon) + (1|exptID), data=onlyCalc)
  summary(testModel8.2)
  rand(testModel8.2)
  r.squaredGLMM(testModel8.2)
  