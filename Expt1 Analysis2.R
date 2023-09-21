#0.open script, set names and parameters
dShell <- read.csv('./shells_Expt1.csv', skip=22)

dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)

#General calculations of mass lost and time
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cTime1 <- strptime(dShell$dateM1, format=dateFormat)
dShell$cTime2 <- strptime(dShell$dateM2, format=dateFormat)
dShell$cDuration <- dShell$cTime2 - dShell$cTime1

#General calculation of shell surface area for each taxa based on generalized shape
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)

dShell$shape <- 'cylinder'
dShell[(dShell$taxon == 'Aragonite'),'shape'] <- 'cube'
dShell[(dShell$taxon == 'Tridacna'),'shape'] <- 'cube'

#pingui and abra technically ovular cylinder
dShell[(dShell$taxon == 'Pinguitellina'), 'shape'] <- 'ovular'
dShell[(dShell$taxon == 'Abranda'), 'shape'] <- 'ovular'

#halimeda different shape entirely... diamond/rhombus
dShell[(dShell$taxon =='Halimeda'), 'shape'] <- 'rhombus'

#basic rectangular cube surface area. does not work for halimeda based off how measurements were taken.
dShell$cSA1 <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#next lines for surface area of cylindrical specimens
subCyl <- which(dShell$shape == 'cylinder')
dShell[subCyl,'cSA1'] <- 2 * (pi* dShell[subCyl,'xDim']/2 * dShell[subCyl,'yDim']/2) + dShell[subCyl,'zDim']*pi*(dShell[subCyl,'yDim'] + dShell[subCyl,'xDim'])/2

#next lines for surface area of ovular specimens
subOv <- which(dShell$shape == 'ovular')
dShell[subOv,'cSA1'] <- 2 * (pi* dShell[subOv,'xDim']/2*dShell[subOv,'yDim']/2) + (pi* dShell[subOv,'zDim']* (dShell[subOv,'xDim']/2 + dShell[subOv,'yDim']/2))

#finally, rough surface area for rhomboid Halimeda - do you need another line to account for the skinny faces?
subRhom <- which(dShell$shape == 'rhombus')
dShell[subRhom,'cSA1'] <- 2 * ((dShell[subRhom,'xDim']*dShell[subRhom,'yDim'])/2)
  
#1.1 just having a look at the measurements to check for errors 

pData <- dShell
TAXA <- unique(pData$taxon)

pdf('./measureCheck.pdf', width=6, height=8, page='A4')

par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(3,3,1,1))

for (t in TAXA) {

	plot(xDim ~ yDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','yDim'))
	plot(xDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'xDim','zDim'))
	plot(yDim ~ zDim, data=pData[(pData$taxon == t),], main=paste(t,'yDim','zDim'))
	plot(xDim ~ thick, data=pData[(pData$taxon == t),], main=paste(t,'xDim','thick'))
	plot(xDim ~ mass1, data=pData[(pData$taxon == t),], main=paste(t,'xDim','mass1'))
	plot(xDim ~ density, data=pData[(pData$taxon == t),], main=paste(t,'xDim','density'))


}

dev.off()

summary(pData)


a <- lm(pMass ~ cSA1  + mass1 + taxon, data=pData)
a <- lm(pMass ~ cSA1  + mass1 + density, data=pData)
summary(a)

a$residuals



#1.2 Total mass lost by taxon separate from variables

pdf('./deltaMass.pdf', width=8, height=6, page='A4')

#edit exptID to be related to this new set of experiments
#for (e in EXPID) {
  
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
#  pData <- dShell[(dShell$exptID == e),]
  #to get anything to work this pData needs to go first!!! added below.
  pData <- dShell
  pData$taxon <- as.factor(pData$taxon)
  taxa <- sort(unique(pData$taxon))
  tFont <- rep(3,length(taxa))
  tFont[which(taxa == 'Aragonite')] <- 1
  #what does following line do? likely have to change integers
  taxaAbrev <- substring(taxa,0,5)
  
  #could skip - just mass lost in mg
  plot(cMass ~ taxon, data=pData[!is.na(pData$cMass),], ylim=c(0,max(pData$cMass, na.rm=TRUE)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData[!is.na(pData$cMass),])
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  #what is the mtext function for
 # mtext(e)


 #1.3 density by taxon
  
  plot(density ~ taxon, data=pData[!is.na(pData$density),], ann=FALSE, axes=FALSE)
  points(density ~ taxon, data=pData[!is.na(pData$density),])
  mtext('Density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)

  
 #1.4 mass lost in percentage
  plot(pMass ~ taxon, data=pData[!is.na(pData$pMass),], ann=FALSE, axes=FALSE)
  points(pMass ~ taxon, data=pData[!is.na(pData$pMass),])
  mtext('Mass lost (%)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
#  mtext(e)
  #not sure what next line is for....
  pData[(pData$pMass>0.3),]
  

#2. Mass lost plotted against surface area
  #2.1 why using cMass (mg) rather than pMass (percentage)????
  if (length(which(!is.na(pData$cSA1))) > 0) {
    plot(cMass/cSA1 ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
    points(cMass/cSA1 ~ taxon, data=pData)
    mtext('Mass lost / SA', side=2, line=3)
    axis(2, las=1)
    axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
 #   mtext(e)
  } else {
    plot(1:1, type='n', ann=FALSE, axes=FALSE)
    mtext('no data', side=1, line=-2)
}

#3. Mass lost plotted against density
  #3.1 Mass/density by taxon
  #this again uses cMass, which is total mass lost, rather than pMass, which is percentage. Change?
if (length(which(!is.na(pData$density))) > 0) {
  plot(cMass/density ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=pData)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
#  mtext(e)
  #following line only needed if you have no data for one part of expt, which you won't
} else {
  plot(1:1, type='n', ann=FALSE, axes=FALSE)
  mtext('no data', side=1, line=-2)
}

dev.off()

#4 Plot these all again - compared to the calculated dimensions and surface area using slightly different calc (flat calc) 

#4.1 Comparison of X and y dimensions by hand vs in imageJ 
pData <- dShell
TAXA <- unique(pData$taxon)

pdf('./measureCheck(measureVScalc).pdf', width=6, height=8, page='A4')
par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(3,3,1,1))

for (t in TAXA) {
plot(xDim ~ calcX, data=pData[(pData$taxon == t),], main=paste(t,'xDim','calcX'))
plot(yDim ~ calcY, data=pData[(pData$taxon == t),], main=paste(t,'yDim','calcY')) }

dev.off()

#4.2 
#general calculation of surface area for each taxa based on flat shape
#do you need to edit next line? 
calcSA <- (pData$calcSA1 + pData$calcSA2) * (1/2) * 2

pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Aragonite')] <- 1
taxaAbrev <- substring(taxa,0,5)
#this flat surface calculation works well for margi, hali, pingui, and abra.
#it should work okay for turbo, arag, tri, and scaph.
#will most likely be worst with liloa, ethalia, and nat.
#but how do we test all of this???

#4.3 first, compare the two calculations of surface area for similarity
pdf('./surfaceAreaCalculation.pdf', width=8, height=6, page='A4')
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))

#can get this to run but it replicates a million times?

plot(cSA1/calcSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='', xlim=c(0,10))
points(cSA1/calcSA ~ taxon, data=pData)
abline(a=1, b=0)
mtext('Caliper SA / Calculated SA', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)


#4.4 second, plot the mass vs calc by taxon

plot(cMass/calcSA ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
points(cMass/calcSA ~ taxon, data=pData)
mtext('Mass lost / Calculated SA', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)

dev.off()

#4.5 then, figuring out which one is better (this or 2.1)! Compare.



#5 Comparison of three replicates for consistency across experiment
  #start by subsetting the three replicates
pData <- dShell
TAXA <- unique(pData$taxon)

pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Aragonite')] <- 1
taxaAbrev <- substring(taxa,0,5)

#this has to come AFTER THE ABOVE
rep1.1 <- subset(pData,exptID =='T1.1')
rep1.2 <- subset(pData,exptID =='T1.2')
rep1.3 <- subset(pData,exptID =='T1.3')

#5.1 Comparison of X, Y, Z, + thickness by taxon per replicate

pdf('./measureCheckReplicates.pdf', width=6, height=8, page='A4')

par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(3,3,1,1))

#how do I plot these but have points from different reps in different shapes?  
  for (t in TAXA) {
    
    plot(xDim ~ yDim, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'xDim','yDim'))
    points(xDim ~ yDim, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(xDim ~ yDim, data=rep1.3[(rep1.3$taxon == t),], pch=17)
    
    plot(xDim ~ zDim, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'xDim','zDim'))
    points(xDim ~ zDim, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(xDim ~ zDim, data=rep1.3[(rep1.3$taxon == t),], pch=17)
    
    plot(yDim ~ zDim, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'yDim','zDim'))
    points(yDim ~ zDim, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(yDim ~ zDim, data=rep1.3[(rep1.3$taxon == t),], pch=17)
    
    plot(xDim ~ thick, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'xDim','thick'))
    points(xDim ~ thick, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(xDim ~ thick, data=rep1.3[(rep1.3$taxon == t),], pch=17)
    
    plot(xDim ~ mass1, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'xDim','mass1'))
    points(xDim ~ mass1, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(xDim ~ mass1, data=rep1.3[(rep1.3$taxon == t),], pch=17)
    
    plot(xDim ~ density, data=rep1.1[(rep1.1$taxon == t),], pch=15, main=paste(t,'xDim','density'))
    points(xDim ~ density, data=rep1.2[(rep1.2$taxon == t),], pch=16)
    points(xDim ~ density, data=rep1.3[(rep1.3$taxon == t),], pch=17)

    }
dev.off()

#5.2 Comparison of SA (caliper) by taxon per replicate

pdf('./ReplicatesMassSA.pdf', width=8, height=6, page='A4')

par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
#for 1.1

  plot(cMass/cSA1 ~ taxon, data=rep1.1, ann=FALSE, axes=FALSE, ylab='', main='Rep1.1')
  points(cMass/cSA1 ~ taxon, data=rep1.1)
  mtext('Mass lost / SA', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)

#for 1.2
  
  plot(cMass/cSA1 ~ taxon, data=rep1.2, ann=FALSE, axes=FALSE, ylab='', main='Rep1.2')
  points(cMass/cSA1 ~ taxon, data=rep1.2)
  mtext('Mass lost / SA', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
#for 1.3
  plot(cMass/cSA1 ~ taxon, data=rep1.3, ann=FALSE, axes=FALSE, ylab='', main='Rep1.3')
  points(cMass/cSA1 ~ taxon, data=rep1.3)
  mtext('Mass lost / SA', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)

  dev.off()  
  
  
#5.3 Comparison of density by taxon per replicate
  pData <- dShell
  TAXA <- unique(pData$taxon)
  
  pdf('./ReplicatesDensity.pdf', width=8, height=6, page='A4')
  
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))

  #attempt to plot them all together but not giving me anything useful
  #SKIP THIS BIT
  plot(cMass/density ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=rep1.1, pch=15)
  points(cMass/density ~ taxon, data=rep1.2, pch=16)
  points(cMass/density ~ taxon, data=rep1.3, pch=17)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #for 1.1
  pdf('./ReplicatesMassDensity.pdf', width=8, height=6, page='A4')
  
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  #for 1.1
  
  plot(cMass/density ~ taxon, data=rep1.1, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=rep1.1)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #for 1.2
  
  plot(cMass/density ~ taxon, data=rep1.2, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=rep1.2)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #for 1.3
  plot(cMass/density ~ taxon, data=rep1.3, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=rep1.3)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  dev.off()  
  
#5.4 Comparison of % mass lost (by taxon?) per replicate
  pdf('./ReplicatesPMassLost.pdf', width=8, height=6, page='A4')
  
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  #for 1.1
  
  plot(pMass ~ taxon, data=rep1.1, ann=FALSE, axes=FALSE, ylab='', ylim=c(0,1))
  points(pMass ~ taxon, data=rep1.1)
  mtext('Mass lost (%)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #for 1.2
  
  plot(pMass ~ taxon, data=rep1.2, ann=FALSE, axes=FALSE, ylab='', ylim=c(0,1))
  points(pMass ~ taxon, data=rep1.2)
  mtext('Mass lost (%)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  #for 1.3
  plot(pMass ~ taxon, data=rep1.3, ann=FALSE, axes=FALSE, ylab='', ylim=c(0,1))
  points(pMass ~ taxon, data=rep1.3)
  mtext('Mass lost (%)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  
  dev.off()  
  
#now altogether???
  
  
#5.5 Compare pH increase over time by replicate (using Hanna logs)
  
pH1.1 <- read.csv('./Expt 1.1 Log 1.csv')
pH1.2 <- read.csv('./Expt 1.2 Log 1.csv')
pH1.3 <- read.csv('./Expt 1.3 Log 1.csv')

#how to combine date and time from two columns into one that makes up time?
dateFormat <- '%m/%d/%y'
timeFormat <- '%H:%M:%S'

#hours since start?
#change to r date (stringtime) and subtract minimum value

#after combining date/time, need to plot date/time ~ pH. shouldn't need ylim if done correctly
plot(Date ~ pH, ylim=c(0,2:32:22), data=pH1.1, type = "o", xlab = "Time", ylab = "pH", main = "Replicate 1.1 pH Over Time")

#repeat for 1.2 and 1.3


#6. Taxon half-life comparison to dissolution rate (only for appropriate taxa)
#work on this later
#subset with Eth, Nat, Turbo, and Pingui