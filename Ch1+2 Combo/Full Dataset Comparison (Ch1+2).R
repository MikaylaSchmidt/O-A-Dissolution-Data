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
dShell1 <- dShell1[,-c(15,17)]
dShell1$chapter <- 'ch1'

#ch2 load
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data")
dShell2 <- read.csv('./Ch2Clean_MasterSet.csv')
dShell2 <- dShell2[,-c(15,16,18,19,20,28, 36,37,38,39,40)]
dShell2$chapter <- 'ch2'

#bind into one
masterSet <- rbind(dShell1, dShell2)


#1.0 t-test to see if similar enough to combine
boxplot(perMassSA~factor(masterSet$chapter)+taxon, data=masterSet)
stripchart(perMassSA ~ factor(masterSet$chapter)+taxon, data = masterSet, vertical=TRUE, add=TRUE, pch=1)

subDouble <- subset(masterSet, masterSet$taxon== 'Aragonite'|masterSet$taxon== 'Calcite' | masterSet$taxon== 'Marginopora')
boxplot(perMassSA~factor(subDouble$chapter)+taxon, data=subDouble, xlab = "Taxon by Expt", ylab='Mass/Surface Area Standardized')
stripchart(perMassSA ~ factor(subDouble$chapter)+taxon, data = subDouble, vertical=TRUE, add=TRUE, pch=1)

#ttests looking at difference here. calc and margi diff is sig
subArag <- subset(subDouble, subDouble$taxon =='Aragonite')
t.test(perMassSA ~ chapter, data=subArag)
subCal <- subset(subDouble, subDouble$taxon== 'Calcite')
t.test(perMassSA ~ chapter, data=subCal)
subMargi <- subset(subDouble, subDouble$taxon== 'Marginopora')
t.test(perMassSA ~ chapter, data= subMargi)


#now do the same but removing ch1expt3 from the mix
subDouble2 <- subset(subDouble, subDouble$exptID != 'T4.1')
subDouble2 <- subset(subDouble2, subDouble2$exptID != 'T4.2')
subDouble2 <- subset(subDouble2, subDouble2$exptID != 'T4.3')

boxplot(perMassSA~factor(subDouble2$chapter)+taxon, data=subDouble2)
stripchart(perMassSA ~ factor(subDouble2$chapter)+taxon, data = subDouble2, vertical=TRUE, add=TRUE, pch=1)

subArag2 <- subset(subDouble2, subDouble2$taxon =='Aragonite')
t.test(perMassSA ~ chapter, data=subArag2)
subCal2 <- subset(subDouble2, subDouble2$taxon== 'Calcite')
t.test(perMassSA ~ chapter, data=subCal2)
subMargi2 <- subset(subDouble2, subDouble2$taxon== 'Marginopora')
t.test(perMassSA ~ chapter, data= subMargi2)


#2.0 looking at data before removing ch1expt3

#looking at all of these in one big data blob
boxplot(perMassSA ~taxon, data=masterSet)

#what are things like if you plot this?
plot(cMass~finalSA, data=masterSet)
plot(percentTotal ~ percentSATotal, data=masterSet, xlab ='Standardized Surface Area', ylab = 'Standardized Mass')
abline(a=0, b=1)
#line of best fit for this
fitLine <- lm(percentTotal ~ percentSATotal, data=masterSet)
abline(fitLine)

#modelling time
model0 <- lm(percentTotal ~ percentSATotal, data=masterSet)
summary(model0)

#what is most important
#THIS DOESN"T RUN
modelFull <- lm(percentTotal ~., data=masterSet)
dredge


#3.0 masterset but with ch1expt3 removed
masterSet2 <- subset(masterSet, masterSet$exptID != 'T4.1')
masterSet2 <- subset(masterSet2, masterSet2$exptID != 'T4.2')
masterSet2 <- subset(masterSet2, masterSet2$exptID != 'T4.3')

#plot to look at the data before modelling
plot(percentTotal ~ percentSATotal, data=masterSet2,xlab ='Standardized Surface Area', ylab = 'Standardized Mass' )
abline(a=0, b=1)
fitLine <- lm(percentTotal ~ percentSATotal, data=masterSet2)
abline(fitLine)

#compare this model to model0
model0.1 <- lm(percentTotal ~ percentSATotal + polymorph, data=masterSet2)
summary(model0.1)

#now begin the process of dredging this like you did with ch 2 in order to figure out 



