dShell <- read.csv('./shells_Expt1.csv', skip=18)

dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)

#General calculations of mass lost and time
dShell$cMass <- -1 * (dShell$mass2 - dShell$mass1)
dShell$pMass <- (dShell$mass1 - dShell$mass2) / dShell$mass1
dShell$cTime1 <- strptime(dShell$dateM1, format=dateFormat)
dShell$cTime2 <- strptime(dShell$dateM2, format=dateFormat)
dShell$cDuration <- dShell$cTime2 - dShell$cTime1

#0.general calculation of shell surface area for each taxa based on generalized shape
#do you need to edit next line? 
dShell$cSize <- (dShell$xDim * dShell$yDim * dShell$zDim) ^ (1/3)

dShell$shape <- 'cylinder'
dShell[(dShell$taxon == 'Aragonite'),'shape'] <- 'cube'
dShell[(dShell$taxon == 'Tridacna'),'shape'] <- 'cube'

#pingui and abra technically ovular cylinder - does that still work?
dShell[(dShell$taxon == 'Pinguitellina'), 'shape'] <- 'ovular'
dShell[(dShell$taxon == 'Abranda'), 'shape'] <- 'ovular'

#halimeda different shape entirely... diamond?

#basic rectangular cube surface area. does not work for halimeda based off how measurements were taken.
dShell$cSA1 <- 2*(dShell$xDim * dShell$yDim) + 2*(dShell$yDim * dShell$zDim) + 2*(dShell$xDim * dShell$zDim)

#next lines for surface area of cylindrical specimens
subCyl <- which(dShell$shape == 'cylinder')
dShell[subCyl,'cSA1'] <- 2 * (pi* dShell[subCyl,'xDim']/2 * dShell[subCyl,'yDim']/2) + dShell[subCyl,'zDim']*pi*(dShell[subCyl,'yDim'] + dShell[subCyl,'xDim'])/2

#next lines for surface area of ovular cylinder specimens
subOv <- which(dShell$shape == 'ovular')
dShell[subOv,'cSA1'] <- 2 * (pi* dShell[subOv,'xDim']*dShell[subOv,'yDim']) + (pi* dShell[subOv,'zDim']* (dShell[subOv,'xDim']* dShell[subOv,'yDim']))
  
#1. Total/Percentage mass lost over time by taxon separate from variables

pdf('./deltaMass.pdf', width=6, height=8, page='A4')

#edit exptID to be related to this new set of experiments
for (e in EXPID) {
  
  par(mfrow=c(3,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  pData <- dShell[(dShell$exptID == e),]
  pData$taxon <- as.factor(pData$taxon)
  taxa <- sort(unique(pData$taxon))
  tFont <- rep(3,length(taxa))
  tFont[which(taxa == 'Aragonite')] <- 1
  #what does following line do? likely have to change integers
  taxaAbrev <- substring(taxa,0,5)
  
  #could skip this following bit? just mass lost in mg
  plot(cMass ~ taxon, data=pData, ylim=c(0,max(pData$cMass)), ann=FALSE, axes=FALSE)
  points(cMass ~ taxon, data=pData)
  mtext('Mass lost (mg)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
  #what is the mtext function for
  mtext(e)
  
  #mass lost in percentage
  plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE)
  points(pMass ~ taxon, data=pData)
  mtext('Mass lost (%)', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  mtext(e)}

#2. Mass lost plotted against surface area
  #why using cMass (mg) rather than pMass (percentage)?
  if (length(which(!is.na(pData$cSA1))) > 0) {
    plot(cMass/cSA1 ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
    points(cMass/cSA1 ~ taxon, data=pData)
    mtext('Mass lost / SA', side=2, line=3)
    axis(2, las=1)
    axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
    mtext(e)
  } else {
    plot(1:1, type='n', ann=FALSE, axes=FALSE)
    mtext('no data', side=1, line=-2)}

#3. Mass lost plotted against density
#this again uses cMass, which is total mass lost, rather than pMass, which is percentage. Change?
if (length(which(!is.na(pData$density))) > 0) {
  plot(cMass/density ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
  points(cMass/density ~ taxon, data=pData)
  mtext('Mass lost / density', side=2, line=3)
  axis(2, las=1)
  axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
  mtext(e)
  #following line only needed if you have no data for one part of expt, which you won't
} else {
  plot(1:1, type='n', ann=FALSE, axes=FALSE)
  mtext('no data', side=1, line=-2)}

#4. Mass lost plotted against both - interaction between two variables
  #can you use the models from the previous two combined with the * interaction? or does that not work for non-linear models?


#Compare above three to fit best model and determine actual significance of each variable
#would it be better to use ANOVA or Chi-squared in this situation? Two indepedent variables, one dependent, categorical...

#5. Plot of pH increase over time for each replicate (not model)
pH1.1 <- read.csv('./1.1 Log 1.csv')

#how to combine date and time?
dateFormat <- '%m/%d/%y'
timeFormat <- '%H:%M:%S'

pH1.1sub <- subset(Date, Time, pH)

plot(pH1.1, type = "o",
     xlab = "Time", ylab = "pH",
     main = "Replicate 1.1 pH Over Time")

#repeat for 1.2 and 1.3

#do you need to use glm for this or just plot the increase over time? not really trying to model anything, just observe...

#6. Taxon half-life comparison to dissolution rate (only for appropriate taxa)
#work on this later
