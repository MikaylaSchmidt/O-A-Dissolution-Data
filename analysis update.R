setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Expt1.csv')
dateFormat <- '%m/%d/%y %H:%M'
EXPID <- unique(dShell$exptID)


pdf('C:/Users/micke/OneDrive/Desktop/deltaMass.pdf', width=6, height=8)

for (e in EXPID) {

	par(mfrow=c(3,1), oma=c(1,1,1,1), mar=c(8,4,1,0))

	pData <- dShell[(dShell$exptID == e),]
	pData$taxon <- as.factor(pData$taxon)
	taxa <- sort(unique(pData$taxon))
	tFont <- rep(3,length(taxa))
	tFont[which(taxa == 'Aragonite')] <- 1
	taxaAbrev <- substring(taxa,0,5)

	plot(cMass ~ taxon, data=pData, ylim=c(0,max(pData$cMass)), ann=FALSE, axes=FALSE)
	points(cMass ~ taxon, data=pData)
	mtext('Mass lost (mg)', side=2, line=3)
	axis(2, las=1)
	axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)
	mtext(e)

	plot(pMass ~ taxon, data=pData, ann=FALSE, axes=FALSE)
	points(pMass ~ taxon, data=pData)
	mtext('Mass lost (%)', side=2, line=3)
	axis(2, las=1)
	axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
	mtext(e)

	if (length(which(!is.na(pData$cSA1))) > 0) {
		plot(cMass/cSA1 ~ taxon, data=pData, ann=FALSE, axes=FALSE, ylab='')
		points(cMass/cSA1 ~ taxon, data=pData)
		mtext('Mass lost / SA', side=2, line=3)
		axis(2, las=1)
		axis(1, at=1:length(taxa), labels=taxa, font=tFont, cex=0.5, las=2)
		mtext(e)
	} else {
		plot(1:1, type='n', ann=FALSE, axes=FALSE)
		mtext('no data', side=1, line=-2)
	}
}

dev.off()

	xVar <- 'xDim'
	yVar <- 'yDim'

#do you really need this bit? what is the point of this?
pdf('./specimenMeasures.pdf', width=6, height=8, page='A4')


for (xVar in c('xDim','yDim','zDim','mass1','SA1'))
for (yVar in c('xDim','mass1','cSize')) {
par(mfrow=c(4,2), oma=c(1,1,1,1), mar=c(8,4,1,0))

for (t in unique(dShell$taxon)) {

	pData <- dShell[(dShell$taxon == t),]
	pData <- pData[!is.na(pData[xVar]),]
	pData <- pData[!is.na(pData[yVar]),]
	
	tFont <- 3
	if (t == 'Aragonite')
		tFont=1

	if (nrow(pData) > 0) {
		plot(pData[,xVar],pData[,yVar], ann=FALSE)
		points(pData[,xVar],pData[,yVar])
		mtext(yVar, side=2, line=3)
		mtext(xVar, side=1, line=3)
		mtext(t, font=tFont)
		a <- lm(pData[,yVar] ~ pData[,xVar])
		rSquared <- round(summary(a)$adj.r.squared,2)
		if (rSquared > 0.6) {
			abline(a)
			mtext(paste('r\u00b2=', rSquared), side=1, adj=0.9, line = -1.5)
		}
	} else {
		plot(1:1, type='n', ann=FALSE, axes=FALSE)
		mtext('no data', side=1, line=-2)
	}
}
}
dev.off()

#for following bit- what does this mean for analysis? why use a t test for this as opposed to others?
#need to work through this bit onwards in order to be able to do anything with it
expMean <- aggregate(dShell$pMass, by=list(dShell$taxon,dShell$exptID),FUN=mean)
expMean <- expMean[(expMean[,'Group.2'] =='P6'),]

expSD <- aggregate(dShell$pMass, by=list(dShell$taxon,dShell$exptID),FUN=sd)
expSD <- expSD[(expSD[,'Group.2'] =='P6'),]

power.t.test(n=5,delta=(expMean[1,'x']-expMean[2,'x']),sd=expSD[1,'x'], sig.level=0.05,power=NULL)


power.t.test(n=5,delta=(expMean[2,'x']-expMean[3,'x']),sd=expSD[3,'x'], sig.level=0.05,power=NULL)
power.t.test(n=NULL,delta=(expMean[7,'x']-expMean[8,'x']),sd=expSD[8,'x'], sig.level=0.05,power=0.8)



dShell <- dShell[(dShell$exptID %in% c('P4','P5','P6')),]

dShell1 <- dShell[(dShell$exptID %in% c('P4')),]

expMean <- aggregate(dShell1$pMass, by=list(dShell1$taxon),FUN=mean)
expSD <- aggregate(dShell1$pMass, by=list(dShell1$taxon),FUN=sd)


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


bob <- data.frame(comparison=tax, power=(res))
bob <- bob[(bob$n > 0),]
bob








