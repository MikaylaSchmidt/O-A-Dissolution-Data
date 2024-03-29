################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2023 Mikayla Schmidt
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
#Fig gen trial

#1.1 seems to actually work to linearize the data for mass1 to pMass relationship
plot(log(mass1)~log(pMass), data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Log of Starting Mass (mg)', xlab='Log of % Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(log(mass1)~log(pMass), data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
#on plotting looks... okayish.
plot(lm.temp)

#abline - separate for non taxon based r2
#r2 is 0.7245 and p value is less than 0.000005
line <- lm(log(mass1)~log(pMass), data=pData)
summary(line)
abline(line)


#1.2 now for finalSA, just logged both again
#halimeda is reallllly weird here
plot(log(finalSA)~log(pMass), data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Surface Area (mm\u00b2)', xlab='% Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(log(finalSA)~log(pMass), data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)
#general abline for trend of taxon independent 
line <- lm(log(finalSA)~log(pMass), data=pData)
summary(line)
abline(line)


#1.3 same but for surface area OVER starting mass
#relationship of surface area to starting mass is a bit confusing... ^(2/3)
plot((finalSA/mass1)~pMass, data=pData, col= pData$tColor, pch=substring(pData$taxon, 0, 2), ylab='Surface Area (mm\u00b2)', xlab='% Mass Lost')
pVal <- vector()
p = 0
for(T in TAXA){
  temp <- pData[(pData$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm((finalSA/mass1)~pMass, data=temp)
  p <- p + 1
  pVal[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',pVal), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#again, but for density
