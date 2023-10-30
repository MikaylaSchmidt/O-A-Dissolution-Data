#Pilot2.1 T Test
dShellFull <- read.csv('./shells_Expt2Pilot.csv', skip=19)
dShell <- subset(dShellFull, exclude == 0)

dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)
dShell$cMassWax <- -1 * (dShell$waxMass2 - dShell$waxMass1)
dShell$waxMassDiff <- dShell$cMassWax-dShell$cMass
dShell$pMassWax <- dShell$cMassWax / dShell$mass1
dShell$wax <- dShell$waxMass1 - dShell$mass1


dShell

plot(cMass~as.factor(taxon), data = dShell)
plot(pMass~as.factor(taxon), data = dShell)
plot(mass1~wax, data=dShell)
plot(cMass~waxMass, data=dShell)
abline(a=0, b=1)
a<- lm(cMass~waxMass, data=dShell)
summary(a)
abline(a)
plot(pMass~pWaxMass, data=dShell)
abline(a=0, b=1)
b<- lm(pMass~pWaxMass, data=dShell)
summary(b)
abline(b)
plot(as.factor(dShell$taxon), dShell$waxMass1-dShell$waxMass2)
points(as.factor(dShell$taxon), dShell$waxMass1-dShell$waxMass2)

plot(cMassWax~as.factor(taxon), data = dShell)
plot(pMassWax~as.factor(taxon), data = dShell)

#final combined pMass and cMass
dShell$cMassFinal <- dShell$cMass
subWax <- which(!is.na(dShell$waxMass))
dShell[subWax, 'cMassFinal'] <- dShell[subWax, 'cMassWax'] 

dShell$pMassFinal <- dShell$pMass 
dShell[subWax, 'pMassFinal'] <- dShell[subWax, 'pMassWax'] 

plot(cMassFinal~as.factor(taxon), data = dShell)
plot(pMassFinal~as.factor(taxon), data = dShell)

#first line runs fine
expMean <- aggregate(dShell$pMassFinal, by=list(dShell$taxon),FUN=mean, na.action=na.omit)
expMean
expMedian <- aggregate(dShell$pMassFinal, by=list(dShell$taxon),FUN=median, na.action=na.omit)
expMedian

expSD <- aggregate(dShell$pMassFinal, by=list(dShell$taxon),FUN=sd)
expSD

power.t.test(n=5,delta=(expMedian[1,'x']-expMedian[2,'x']),sd=expSD[2,'x'], sig.level=0.05,power=NULL)


power.t.test(n=5,delta=(expMean[2,'x']-expMean[3,'x']),sd=expSD[3,'x'], sig.level=0.05,power=NULL)
power.t.test(n=NULL,delta=(expMedian[5,'x']-expMedian[4,'x']),sd=expSD[4,'x'], sig.level=0.05,power=0.8)




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
      
      #		res[c] <- power.t.test(n=NULL,delta=DELTA,sd=SD, sig.level=0.05,power=0.8)$n
      res[c] <- power.t.test(n=5,delta=DELTA,sd=SD, sig.level=0.05,power=NULL)$power
    }
  }}

bob <- data.frame(comparison=tax, n=round(res))
bob <- bob[(bob$n > 0),]
bob
