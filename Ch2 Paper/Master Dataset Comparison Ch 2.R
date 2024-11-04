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
  
  #fixing issues with factoring
  TAXA <- unique(masterSet$taxon)
  masterSet$taxon <- as.factor(masterSet$taxon)
  taxa <- sort(unique(masterSet$taxon))
  tFont <- rep(3,length(taxa))
  tFont[which(taxa == 'Anadara')] <- 1
  taxaAbrev <- substring(taxa,0,5)
  TAXA2 <- data.frame(taxon=taxa, tFont = tFont, tPoint=substring(taxa,0,1), tAbrev = taxaAbrev, tColor=rainbow(length(taxa)), tLine = 1:length(taxa))
  masterSet<-merge(masterSet,TAXA2,by='taxon')
  
  
  #porosity  calc
  masterSet$poros <- masterSet$mass1/(masterSet$volume * 2.83)
  subCal <- which(masterSet$polymorph =='Calcite')
  masterSet[subCal, 'poros'] <- masterSet[subCal, 'mass1'] / (masterSet[subCal, 'volume'] * 2.711)
  masterSet$poros <- 1- masterSet$poros
  
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


#quick save for later
write.csv(masterSet, 'Ch2Clean_MasterSet.csv', row.names=TRUE)
dev.off()
  
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
  #masterSet <- subset(masterSet, exptNo != 'Expt2')
  
  plot(massSA~taxon, data=masterSet)
  points(massSA~taxon, data=masterSet)
  factor(masterSet$exptNo)
  boxplot(massSA~factor(masterSet$exptNo)+taxon, data=masterSet)
  stripchart(massSA ~ factor(masterSet$exptNo)+droplevels(taxon), data = masterSet, vertical=TRUE, add=TRUE, pch=1)
  str(masterSet)
  
  plot(perMassSA~taxon, data=masterSet)
  points(perMassSA~taxon, data=masterSet)
  boxplot(perMassSA~factor(masterSet$exptNo)+taxon, data=masterSet)
  stripchart(perMassSA ~ factor(masterSet$exptNo)+droplevels(taxon), data = masterSet, vertical=TRUE, add=TRUE, pch=1)
  
  
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
  summary1_3$n <- 14
  #sub14 <- which(summary1_3$taxon == 'Marginopora'|summary1_3$taxon=='Turbo')
  #summary1_3[sub14, 'n'] <-14
  #sub13 <- which(summary1_3$taxon == 'Alaona')
  #summary1_3[sub13, 'n'] <-13
  #sub12 <- which(summary1_3$taxon == 'Liloa')
  #summary1_3[sub12, 'n'] <-12
  
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
  
  
#3.0 Now onto %/% standard
  
  #set colours using this line here
  myCol <- ifelse(levels(factor(Expt1_3$exptNo))=="Expt1" , rgb(0.1,0.1,0.7,0.5) , 
                  ifelse(levels(factor(Expt1_3$exptNo))=="Expt3", rgb(0.8,0.1,0.3,0.6),
                         "grey90" ) )
  
  #boxplot of the things you actually want to look at
  boxplot(perMassSA~exptNo * droplevels(taxon), data= Expt1_3, ylab = 'Mass/Surface Area standardized', xlab = 'Taxon by Expt', col= myCol, xaxt ='n')
  axis(1, at = c(1.5,3.5,5.5,7.5,9.5,11.5), labels=taxa2, las=1)
  stripchart(perMassSA ~ factor(Expt1_3$exptNo)+droplevels(taxon), data = Expt1_3, vertical=TRUE, add=TRUE, pch=1)
  abline(h=1)
  
  #next bit is trial that requires percentTotal/percent SA total calc  
  Expt1_3$perMassSA <- Expt1_3$percentTotal/Expt1_3$percentSATotal
  summary1_3Per <- data.frame(taxon = TAXA2, p.value = NA, a.p.value = NA, mean.diff=NA, conf.int1 = NA, conf.int2 = NA , n = 17, SD = NA, power = NA)
  sub15 <- which(summary1_3Per$taxon == 'Saccostrea')
  summary1_3Per[sub15, 'n'] <-15
  sub7 <- which(summary1_3Per$taxon == 'Centrostephanus')
  summary1_3Per[sub7, 'n'] <-7
  
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
    summary1_3Per[p, 'SD'] <- round(ttest.temp$p.value, digits=4)
    summary1_3Per[p, 'p.value'] <- round(ttest.temp$p.value, digits=4)
    summary1_3Per[p, 'power'] <- round(power.temp$power, digits=2)
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
  plot(perMassSA~taxon, data=masterSet, xlab = 'Taxon', ylab = 'Standardized Mass Lost/Surface Area')
  points(perMassSA~taxon, data=masterSet)
  abline(h=1)

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
        ttest.temp <- t.test(perMassSA ~ taxon, data = temp3)
        res[c] <- ttest.temp$p.value
      }
    }}
  
  dataTrial <- data.frame(comparison=tax, p.value=round(res, digits=6), adj.p.value =NA)
  dataTrial <- dataTrial[(dataTrial$p.value != -99),]
  dataTrial$adj.p.value <- p.adjust(dataTrial$p.value)
  dataTrial
  
  
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
  
  
  
#6.0 Linear modelling time!
  #a few linear models to compare unstandardized vs standardized with the same measurements
  
  #non-standardized
    plot(cMass~finalSA, data= Expt123)
    model0 <- lm(cMass~finalSA, data= Expt123)
    summary(model0)
    model1 <- lm(cMass~finalSA + taxon, data= Expt123)
    summary(model1)
    model2 <- lm(cMass~finalSA + polymorph, data= Expt123)
    summary(model2)
  
  #better r2 for standardized
    model3 <- lm(percentTotal~percentSATotal, data= Expt123)
    summary(model3)
    model4 <- lm(percentTotal~percentSATotal + taxon, data= Expt123)
    summary(model4)

  
#7.0 On to model selection as a whole, no taxon version for model selection
#create new data frame using only the relevant data aka relData
relData <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptNo', 'exptID', 'deviationStan', 'taxon', 'polymorph', 'volume')]
noTaxon <- relData[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptNo', 'exptID', 'deviationStan', 'polymorph', 'volume')]

#dredge function load in
library(lmerTest)
library(MuMIn)

#7.1 first model: dredging up best model without taxon
  model1 <- lm(percentTotal ~., data= noTaxon)
  options(na.action = 'na.fail')
  dredge(model1)
  summary(model1)
  AICc(model1)
  
#the best model without taxon, as dredged previously (doesn't deal with the exptID issue)
#also exptID is significant in dredge, but not sig in modelBest1 & doing nothing for r2 so should be taken out
  modelBest1 <- lm(percentTotal ~ percentSATotal + cSize + polymorph + exptID, data= noTaxon)
  summary(modelBest1)
  AICc(modelBest1)
  modelBest1.2 <- lm(percentTotal ~ percentSATotal + cSize + polymorph, data= noTaxon)
  summary(modelBest1.2)
  AICc(modelBest1.2)

#7.2 second model: dredging up beset model with taxon
  model0 <- lm(percentTotal ~., data= relData)
  summary(model0)
  AICc(model0)
  dredge(model0, evaluate = TRUE)
  summary(model0)
#the best model with taxon, as dredged, only includes taxon and %total
  modelBest0 <- lm(percentTotal ~ percentSATotal + taxon, data= relData)
  dredge(modelBest0)
  summary(modelBest0) 
  AICc(modelBest0)

  
#8.0 Doing the same thing as previously but using mixed linear model instead of dredge fxn
library(cAIC4)
model0.1 <- lmer(percentTotal~percentSATotal + cSize + poros + deviationStan  + (1|exptNo) + (1|taxon), data = relData)
dredge(model0.1)

  #looking at same model with taxon vs with polymorph and potentially mag
  #interestingly this seems to indicate polymorph and taxon have roughly the same strength and cSize itself is sig
  summary(lm(percentTotal~percentSATotal, data = relData))
  r.squaredGLMM(lmer(percentTotal~percentSATotal + (1|taxon), data = relData))
  r.squaredGLMM(lmer(percentTotal~percentSATotal + (1|polymorph), data = relData))
  r.squaredGLMM(lmer(percentTotal~percentSATotal + (1|polymorph) + (1|taxon), data = relData))

#quick try of a stepwise regression just for the fun of it  
stepmod1 <- lm(percentTotal ~ ., data=relData)
  step(stepmod1, direction =c('both'))
stepmod2 <- lm(percentTotal ~ ., data=noTaxon)
  step(stepmod2, direction =c('both'))
summary(lm(percentTotal ~ percentSATotal + cSize + exptNo + 
       deviationStan + polymorph, data = noTaxon))


#7.0 Comparing stuff - what if you look at just calcite specs? Does that improve things???????
#first, calcite subset
onlyCalc <- subset(masterSet, masterSet$polymorph == 'Calcite')
droplevels(as.factor(onlyCalc$polymorph))
droplevels(as.factor(onlyCalc$taxon))
onlyCalc$taxon <- factor(onlyCalc$taxon, levels = c('Calcite',  'Centrostephanus', 'Pecten', 'Saccostrea', 'Marginopora'))


#plot of massSA
boxplot(massSA ~ taxon, data = onlyCalc)
points(massSA~taxon, data=onlyCalc)

#is mass/SA plot better fit by only calcite? no.ignore this bit in chapter.
plot(cMass ~ finalSA, data=onlyCalc)
only1 <- lm(cMass ~ finalSA, data=onlyCalc)
summary(only1)
only2 <- lm(cMass ~ finalSA * mag, data=onlyCalc)
summary(only2)

#plot of perMassSA
boxplot(perMassSA ~taxon, data = onlyCalc, xlab='Taxon', ylab ="Mass/Surface Area (Standardized)")
points(perMassSA~taxon, data=onlyCalc)

#perMassSA plot fit with just calcite 
#this is somehow not actually that much better than the previous
#probably because saccostrea is shit
plot(percentTotal ~ percentSATotal, data=onlyCalc, ylab = 'Standardized Mass', xlab='Standardized Surface Area')
only3 <- lm(percentTotal ~ percentSATotal, data=onlyCalc)
summary(only3)
only4 <- lm(percentTotal ~ percentSATotal + taxon, data=onlyCalc)
summary(only4)


#no saccostrea? this is actually really very good r2. have more look at the plot later
noSacc <- subset(onlyCalc, onlyCalc$taxon != 'Saccostrea')
boxplot(massSA ~taxon, data = noSacc)
plot(percentTotal ~ percentSATotal, data=noSacc)
only5 <- lm(percentTotal ~ percentSATotal, data=noSacc)
summary(only5)
only6 <- lm(percentTotal ~ percentSATotal * mag, data=noSacc)
summary(only6)
