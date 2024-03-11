#Ch 1 Fig Gen 2 trial script
#for things you don't want but should keep

#3.2.2 shell mean dissolution by taxon
Expt1_3$meanDiss <- ((Expt1_3$cMass)^(2/3))/48 
sub3 <- which(Expt1_3$exptNo == 'Expt3')
Expt1_3[sub3,'meanDiss'] <- Expt1_3[sub3,'cMass']^(2/3)/168

plot(meanDiss ~ taxon, data=Expt1_3, ann=FALSE, axes=FALSE)
points(meanDiss ~ taxon, data=Expt1_3)
mtext('Mean Dissolution per Hour', side=2, line=3)
axis(2, las=1)
axis(1, at=1:length(taxa), labels=taxa, font=3, cex=0.5, las=2)

plot(meanDiss~mass1, data=Expt1_3, col= Expt1_3$tColor , pch=substring(Expt1_3$taxon, 0, 2))
r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA){
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(meanDiss~mass1, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)


#Analysis to examine % dissolutions by % volume (of total)
#again but for volume instead of SA 
volumeTotal <- aggregate(Expt1_3$volume, by = list(Expt1_3$exptID),FUN=sum)
colnames(volumeTotal)<-c('exptID','volumeTotal')
Expt1_3 <- merge(Expt1_3, volumeTotal, by = 'exptID')
Expt1_3$percentVolumeTotal <- (Expt1_3$volume)/(Expt1_3$volumeTotal) * 100

plot(percentTotal~percentVolumeTotal, data = Expt1_3, col= Expt1_3$tColor, pch=substring(Expt1_3$taxon, 0, 2), xlab = '% of Total Volume', ylab = "% of Total Mass Lost")
line <-lm(percentTotal~percentVolumeTotal, data = Expt1_3)
abline(lm(percentTotal~percentVolumeTotal, data = Expt1_3))
abline(a=0, b=1)

r2 <- vector()
slope <- vector()
p = 0
for(T in TAXA2$taxon){
  temp <- Expt1_3[(Expt1_3$taxon == T),]
  pref <- TAXA2[(TAXA2$taxon==T),]
  lm.temp <- lm(percentTotal~percentVolumeTotal, data=temp)
  p <- p + 1
  r2[p] <- round(summary(lm.temp)$adj.r.squared,3)
  slope[p] <- round(coefficients(lm.temp)[2],3)
  abline(lm.temp, col=pref$tColor, lty=pref$tLine)
}
legend('topright', legend= paste(TAXA2$taxon,' r\u00b2=',r2), col=TAXA2$tColor, lty=TAXA2$tLine, text.font = TAXA2$tFont, pch= TAXA2$tPoint, text.col= TAXA2$tColor)

dev.off()

#SA seems to be 1:1 while volume is closer to 1:2
Expt1_3$newRatio <- (Expt1_3$percentTotal)/(Expt1_3$percentVolumeTotal)
plot(newRatio~taxon, data=Expt1_3, ylab='% Mass Lost/% Volume')
abline (a=1, b=0)
