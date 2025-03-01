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
#Old Script Full Ch2 Model
#Read in first half of 'Master Dataset Comparison Ch 2' first, up to end of part 5.0

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
relData <- Expt123[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptID', 'deviationStan', 'taxon', 'polymorph', 'volume')]
noTaxon <- relData[,c('percentTotal', 'percentSATotal', 'cSize', 'poros', 'exptID', 'deviationStan', 'polymorph', 'volume')]

#dredge function load in
library(lmerTest)
library(MuMIn)

#visual of the plot
plot(percentTotal ~ percentSATotal, data= relData)
fitModel <- lm(percentTotal ~ percentSATotal, data= relData)
plot(fitModel)

#based off of this, I think you are going to have to log the data
plot(log(percentTotal) ~ log(percentSATotal), data= relData)
modelLog <- lm(log(percentTotal) ~ log(percentSATotal), data= relData)
plot(modelLog)
summary(modelLog)

modelLog2 <- lm(log(cMass)~ log(finalSA), data=Expt123)
summary(modelLog2)


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

#6.0 Doing the same thing as previously but using mixed linear model instead of dredge fxn

model0.1 <- lmer(percentTotal~percentSATotal + cSize + poros + deviationStan  + (1|exptNo) + (1|taxon), data = relData)
dredge(model0.1)

#looking at same model with taxon vs with polymorph and potentially mag
#interestingly this seems to indicate polymorph and taxon have roughly the same strength and cSize itself is sig
summary(lm(percentTotal~percentSATotal, data = relData))
summary(lm(log(percentTotal)~log(percentSATotal), data = relData))
r.squaredGLMM(lmer(percentTotal~percentSATotal + (1|exptID), data = relData))
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

