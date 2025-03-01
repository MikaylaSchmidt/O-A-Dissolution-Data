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
#Ch 1 + 2 Master Dataset Comparison
#for looking at allllll the data together

#ch1 load
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell1 <- read.csv('./Ch1Clean_MasterSet.csv')
dShell1 <- dShell1[,-c(1,15,23)]
dShell1$chapter <- 'ch1'


#ch2 load
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell2 <- read.csv('./Ch2Clean_MasterSet.csv')
dShell2 <- dShell2[,-c(1, 15,16,19,20,23,24,25,28,32)]
dShell2$chapter <- 'ch2'
dShell2$waxYN <- 'No Wax'

#bind into one
masterSet <- rbind(dShell1, dShell2)


#1.0 t-test to see if similar enough to combine
par(mfrow=c(1,1))
boxplot(perMassSA~factor(masterSet$chapter)+taxon, data=masterSet)
stripchart(perMassSA ~ factor(masterSet$chapter)+taxon, data = masterSet, vertical=TRUE, add=TRUE, pch=1)

subDouble <- subset(masterSet, masterSet$taxon== 'Aragonite'|masterSet$taxon== 'Calcite' | masterSet$taxon== 'Marginopora')
boxplot(perMassSA~factor(subDouble$chapter)+taxon, data=subDouble, xlab = "Taxon by Chapter", ylab='Mass/Surface Area Standardized')
stripchart(perMassSA ~ factor(subDouble$chapter)+taxon, data = subDouble, vertical=TRUE, add=TRUE, pch=1)

boxplot(massSA~factor(subDouble$chapter)+taxon, data=subDouble, xlab = "Taxon by Chapter", ylab='Mass/Surface Area')
stripchart(massSA ~ factor(subDouble$chapter)+taxon, data = subDouble, vertical=TRUE, add=TRUE, pch=1)


#ttests looking at difference here. calc and margi diff is sig, this is why you need to remove expt3
#just using massSA but can use perMassSA but have to p-adjust that
Ch1_2 <- data.frame(taxon= c('Aragonite', 'Calcite', 'Marginopora'), p.value =NA, adj.p.value =NA)

subArag <- subset(subDouble, subDouble$taxon =='Aragonite')
shapiro.test(subArag$massSA)
t.test(massSA ~ chapter, data=subArag)
wilcox.test(massSA ~ chapter, data=subArag)
Ch1_2[1,2]<- wilcox.test(massSA ~ chapter, data=subArag)$p.value

subCal <- subset(subDouble, subDouble$taxon== 'Calcite')
t.test(massSA ~ chapter, data=subCal)
wilcox.test(massSA ~ chapter, data=subCal)
Ch1_2[2,2] <- wilcox.test(massSA ~ chapter, data=subCal)$p.value

subMargi <- subset(subDouble, subDouble$taxon== 'Marginopora')
t.test(massSA ~ chapter, data= subMargi)
wilcox.test(massSA ~ chapter, data=subMargi)
Ch1_2[3,2] <- wilcox.test(massSA ~ chapter, data=subMargi)$p.value

Ch1_2$adj.p.value <- p.adjust(Ch1_2$p.value)
print(Ch1_2)
xtable(Ch1_2)
#now do the same but removing ch1expt3 from the mix
subDouble2 <- subset(subDouble, subDouble$exptID != 'T3.1')
subDouble2 <- subset(subDouble2, subDouble2$exptID != 'T3.2')
subDouble2 <- subset(subDouble2, subDouble2$exptID != 'T3.3')

boxplot(massSA~factor(subDouble2$chapter)+taxon, data=subDouble2, xlab= 'Taxon', ylab='Mass Lost/Surface Area', col=c('#CC1A4D99', '#1A1AB380'), xaxt='n')
stripchart(massSA ~ factor(subDouble2$chapter)+taxon, data = subDouble2, vertical=TRUE, add=TRUE, pch=1)
mtext('Aragonite', side=1, line=1, at=1.5)
mtext('Calcite', side=1, line=1, at=3.5)
mtext(~italic('Marginopora'), side=1, line=1.1, at=5.5)
legend(x= 'topleft', legend = c('Chapter 1','Chapter 2'), col= c('#CC1A4D99', '#1A1AB380'), cex=1.2, pch = 15)


subArag2 <- subset(subDouble2, subDouble2$taxon =='Aragonite')
shapiro.test(subArag2$perMassSA)
#t.test(perMassSA ~ chapter, data=subArag2)
wilcox.test(massSA ~ chapter, data=subArag2)
Ch1_2[1,2]wilcox.test(perMassSA ~ chapter, data=subArag2)

subCal2 <- subset(subDouble2, subDouble2$taxon== 'Calcite')
shapiro.test(subCal2$perMassSA)
#t.test(perMassSA ~ chapter, data=subCal2)
wilcox.test(massSA ~ chapter, data=subCal2)
wilcox.test(perMassSA ~ chapter, data=subCal2)

subMargi2 <- subset(subDouble2, subDouble2$taxon== 'Marginopora')
shapiro.test(subMargi2$perMassSA)
#t.test(perMassSA ~ chapter, data= subMargi2)
wilcox.test(massSA ~ chapter, data=subMargi2)
wilcox.test(perMassSA ~ chapter, data=subMargi2)


#2.0 looking at data before removing ch1expt3
#skip all of this for final chapter

#looking at all of these in one big data blob
boxplot(perMassSA ~taxon, data=masterSet)

#what are things like if you plot this?
#if you log this data, the ch1 expt3 margi and calcite become very clearly weird
plot(cMass~finalSA, data=masterSet)
plot(log(cMass)~log(finalSA), data=masterSet)
plot(percentTotal ~ percentSATotal, data=masterSet, xlab ='Standardized Surface Area', ylab = 'Standardized Mass')
abline(a=0, b=1)
#line of best fit for this
fitLine <- lm(percentTotal ~ percentSATotal, data=masterSet)
abline(fitLine)

#modelling time
model0 <- lm(percentTotal ~ percentSATotal, data=masterSet)
summary(model0)



#3.0 masterset but with ch1expt3 removed
masterSet2 <- subset(masterSet, masterSet$exptID != 'T3.1')
masterSet2 <- subset(masterSet2, masterSet2$exptID != 'T3.2')
masterSet2 <- subset(masterSet2, masterSet2$exptID != 'T3.3')

#plot to look at the data before modelling
plot(cMass ~ finalSA, data = masterSet2)
plot(log(cMass) ~ log(finalSA), data = masterSet2)
plot(percentTotal ~ percentSATotal, data=masterSet2,xlab ='Standardized Surface Area', ylab = 'Standardized Mass' )
abline(a=0, b=1)
fitLine <- lm(percentTotal ~ percentSATotal, data=masterSet2)
abline(fitLine)

#compare this model to model0
#most likely going to need to log this data like previous
model0.1 <- lm(percentTotal ~ percentSATotal, data=masterSet2)
summary(model0.1)
plot(model0.1)

#second model but logged
model0.2 <- lm(log(percentTotal) ~ log(percentSATotal), data=masterSet2)
summary(model0.2)
plot(model0.2)


#3.0 now begin the process of dredging this like you did with ch 2 in order to figure out 
library(lmerTest)
library(MuMIn)

relDataFull <- masterSet2[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptID', 'deviationStan', 'taxon', 'polymorph')]
noTaxonFull <- relDataFull[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptID', 'deviationStan', 'polymorph')]
options(na.action = 'na.fail')

#first dredging the taxon full model, compare to model0.1
model1 <- lm(percentTotal ~ ., data=relDataFull)
summary(model1)
dredge(model1)

model2 <- lm(percentTotal ~ percentSATotal + poros + taxon, data=relDataFull)
summary(model2)

#second dredging the taxon full model, compare to model0.1 and model1
model3 <- lm(percentTotal ~ ., data=noTaxonFull)
summary(model3)
dredge(model3)


#4.0 making visuals of the model data
