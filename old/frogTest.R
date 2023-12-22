#frog test
library(AICcmodavg)
data("dry.frog")

##extract only first 7 columns
frog <- dry.frog[, 1:7]

##first lines
head(frog)
str(frog)
any(is.na(frog))

#1.Null
#2.Shade
#3.Substrate
#4.Shade + substrate
#5.Null and mass
#6.Shade and mass
#7.Subrate and mass
#8.Shade + substrate with mass

##center initial mass
#i don't quite understand the purpose of doing this
frog$InitMass_cent <- frog$Initial_mass - mean(frog$Initial_mass)
frog$InitMass2 <- frog$InitMass_cent^2

##run global model
global <- lm(Mass_lost ~ InitMass_cent + InitMass2 + Substrate + Shade,
               data = frog)
par(mfrow = c(2, 2))
plot(global)

frog$logMass_lost <- log(frog$Mass_lost + 1) 
#adding 1 due to presence of 0's

#log global model
global.log <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Substrate + Shade,
                 data = frog)
par(mfrow = c(2, 2))
plot(global.log)

#fit all the models
m.null <- lm(logMass_lost ~ 1,
             data = frog)
m.shade <- lm(logMass_lost ~ Shade,
                data = frog)
m.substrate <- lm(logMass_lost ~ Substrate,
                    data = frog)
m.shade.substrate <- lm(logMass_lost ~ Shade + Substrate,
                          data = frog)
m.null.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2,
                    data = frog)
m.shade.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Shade,
                     data = frog)
m.substrate.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Substrate,
                         data = frog)
m.global.mass <- global.log

##store models in named list
Cand.models <- list("null" = m.null, "shade" = m.shade,
                      "substrate" = m.substrate,
                      "shade + substrate" = m.shade.substrate,
                      "mass" = m.null.mass, "mass + shade" = m.shade.mass,
                      "mass + substrate" = m.substrate.mass,
                      "global" = m.global.mass)
selectionTable <- aictab(cand.set = Cand.models)
selectionTable

##can export table in LaTeX for formatting reasons

##confidence set of models to summarize what is best 
confset(cand.set = Cand.models)

#then confidence intervals
confint(m.global.mass)


modavg(cand.set = Cand.models, parm = "Shade")
modavg(Cand.models, parm = "SubstrateSPHAGNUM")
modavgShrink(cand.set = Cand.models, parm = "Shade")

