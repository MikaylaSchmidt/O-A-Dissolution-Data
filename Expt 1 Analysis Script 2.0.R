#New Script Data
library(AICcmodavg)
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')

x <- dShell[c(9,25,26,36,37,38)]
cor(log(x[,c(1,4,5)]))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1)))

summary(lm(log(x$cMass)~log(x$mass1)+log(x$finalSA)+log(x$volume)))
summary(lm(log(x$cMass)~log(x$mass1)+log(x$finalSA)))
summary(lm(log(x$cMass)~log(x$mass1)+log(x$volume)))
summary(lm(log(x$cMass)~log(x$mass1)))

summary(lm(log(x$cMass)~log(x$mass1/x$volume)))
summary(lm(log(x$mass1-x$cMass)~log(x$mass1/x$volume)))

x2 <- dShell[c(9,16,25,26,36,37,38)]
summary(lm(log(x2$mass1-x2$cMass)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))
summary(lm(log(x2$cMass)~log(x2$mass1)+log(x2$finalSA)+log(x2$volume)+log(x2$thick)))

minilso(log(x2$cMass),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))
minilso(log(x2$mass1 - x2$cMass),log(cbind(x2$mass1, x2$thick, x2$finalSA, x2$volume)))

x3 <- data.frame(scale(log(dShell[c(9,16,25,26,36,37,38)])))
summary(lm((x3$cMass)~(x3$mass1)+(x3$finalSA)+(x3$volume)+(x3$thick)))
x3

#mixed model

#linear mixed effects -lmer
#lmerTest
library(lmerTest)
summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))

summary(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
rand(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
ranova(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
 
#but no R squared until   
library(MuMIn)
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$cMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

#now time for the Percent Mass
hist(dShell$pMass)  
hist(log(1-dShell$pMass))

r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$exptID)))
r.squaredGLMM(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon)))

summary(lmer(x3$pMass~ x3$mass1 + x3$finalSA + x3$volume + x3$thick + (1|dShell$taxon) + (1|dShell$exptID)))
