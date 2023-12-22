#frogTestShell
library(AICcmodavg)
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')

#is this next bit important/relevant
sumTable <- aggregate(dShell$taxon, by=list(dShell$taxon), FUN=length)
sumTable
colnames(sumTable) <- c('taxon', 'n')
sumTable$pValue = 'na'

#full data model
fullDataModel <- dShell[c(3,4,9,25,26,36,37,38)]
fullDataModel

#1.1 This is for cMass.

#split into dataset with important variables
cData <- dShell[c(3,4,9,25,36,37,38)]
head(cData)
str(cData)
any(is.na(cData))

#0. Global hypothesis
#Hypothesis: Mass lost by shells varies with surface area, volume, density, and mass.

#first, need to center initial mass (aka variable mass1)
cData$InitMass_cent <- cData$mass1 - mean(cData$mass1)
#then square it
cData$InitMass2 <- cData$InitMass_cent^2

#then, the global model itself.....
global <- lm(cMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
             data = cData)
par(mfrow = c(2, 2))
plot(global)

#try this next - log transform mass lost
cData$log_cMass <- log(cData$cMass + 1)

#run global model again
global.log <- lm(log_cMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
                 data = cData)
par(mfrow = c(2, 2))
plot(global.log)

#1. Null
#Hypothesis: Mass lost by shells is constant.
m.null <- lm(log_cMass ~ 1, data = cData)

#2. finalSA
#Hypothesis: Mass lost by shells varies with surface area.
m.finalSA <- lm(log_cMass ~ finalSA, data = cData)

#3. volume
#Hypothesis: Mass lost by shells varies with volume.
m.volume <- lm(log_cMass ~ volume, data = cData)

#4. densityMV
#Hypothesis: Mass lost by shells varies with density.
m.densityMV <- lm(log_cMass ~ densityMV, data = cData)

#5. finalSA + volume
#Hypothesis: Mass lost by shells varies with surface area and volume.
m.finalSA.volume <- lm(log_cMass ~ finalSA + volume, data = cData)

#6. finalSA + densityMV
#Hypothesis: Mass lost by shells varies with surface area and density.
m.finalSA.densityMV <- lm(log_cMass ~ finalSA + densityMV, data = cData)

#7. volume + densityMV
#Hypothesis: Mass lost by shells varies with volume and density.
m.volume.densityMV <- lm(log_cMass ~ volume + densityMV, data = cData)

#8. Null and mass
#Hypothesis: Mass lost by shells varies with shell initial mass.
m.null.mass <- lm(log_cMass ~ InitMass_cent + InitMass2, data = cData)

#9. finalSA and mass
#Hypothesis: Mass lost by shells varies with shell surface area and initial mass. 
m.finalSA.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + finalSA, data = cData)

#10. volume and mass
#Hypothesis: Mass lost by shells varies with shell volume and initial mass.
m.volume.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + volume, data = cData)

#11.density and mass
#Hypothesis: Mass lost by shells varies with shell density and initial mass.
m.densityMV.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + densityMV, data = cData)

#12. finalSA + volume and mass
#Hypothesis: Mass lost by shells varies with surface area, volume and mass.
m.finalSA.volume.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + finalSA + volume, data = cData)

#13. finalSA + densityMV and mass
#Hypothesis: Mass lost by shells varies with surface area, density, and mass.
m.finalSA.densityMV.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + finalSA + densityMV, data = cData)

#14. volume + densityMV and mass
#Hypothesis: Mass lost by shells varies with volume, density, and mass. 
m.volume.densityMV.mass <- lm(log_cMass ~ InitMass_cent + InitMass2 + volume + densityMV, data = cData)

#15. Global Model
m.global.mass <- global.log

##store models in named list
Cand.models <- list("null" = m.null, "surface area" = m.finalSA,
                      "volume" = m.volume, "density" = m.densityMV,
                      "surface area + volume" = m.finalSA.volume,
                      "surface area + density" = m.finalSA.densityMV,
                      "volume + density" = m.volume.densityMV,
                      "mass" = m.null.mass, 
                      "mass + surface area" = m.finalSA.mass,
                      "mass + volume" = m.volume.mass,
                      "mass + density" = m.densityMV.mass,
                      "mass + surface area + volume" = m.finalSA.volume.mass,
                      "mass + surface area + density" = m.finalSA.densityMV.mass,
                      "mass + volume + density" = m.volume.densityMV.mass,
                      "global" = m.global.mass)

#model selection table time! using AICc rather than AIC
selectionTable <- aictab(cand.set = Cand.models)
selectionTable

#nice pretty table time
library(xtable)
print(xtable(selectionTable, caption = "Model selection table on shell mass lost.",
               label = "tab:selection"),
        include.rownames = FALSE, caption.placement = "top")

#2.1 Now..... do it all again for pMass!

#split into dataset with important variables
pData <- dShell[c(3,4,9,26,36,37,38)]
head(pData)
str(pData)
any(is.na(pData))

#0. Global hypothesis
#Hypothesis: Mass lost by shells varies with surface area, volume, density, and mass.

#first, need to center initial mass (aka variable mass1)
pData$InitMass_cent <- cData$mass1 - mean(cData$mass1)
#then square it
pData$InitMass2 <- cData$InitMass_cent^2

#then, the global model itself.....
global <- lm(pMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
             data = pData)
par(mfrow = c(2, 2))
plot(global)

#try this next - log transform mass lost
pData$log_pMass <- log(pData$pMass + 1)

#run global model again
global.log <- lm(log_pMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
                 data = pData)
par(mfrow = c(2, 2))
plot(global.log)

#ok that's actually still crap on the Q-Q. how about sq root?
pData$sqrt_pMass <- sqrt(pData$pMass)
global.sqrt <- lm(sqrt_pMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
                 data = pData)
par(mfrow = c(2, 2))
plot(global.sqrt)

#how about..... cube root?
pData$cuberoot_pMass <- (pData$pMass)^(1/3)
global.cube <- lm(cuberoot_pMass ~ InitMass_cent + InitMass2 + finalSA + volume + densityMV,
                  data = pData)
par(mfrow = c(2, 2))
plot(global.cube)

#that is literally not better at all, so come back to this on monday.
