#ggplot dissolution paste
#removed from 'Ch 1 Fig Generation 2 on 3/26/24

#3.2 Setup and execution of comparison plots
#only thing necessary for this is to split it into two more subsets by overall experiment 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 Paper Figs")
pdf('Expt1&3Comp.pdf', height = 6 , width = 8)
par(mfrow=c(1,1), oma=c(1,1,1,1))

#firstly, some setup to subset
library(ggplot2)
factor(Expt123$exptNo)
Expt123$exptNo <- as.factor(Expt123$exptNo)
str(Expt123)

#cMass by taxon, split by experiment
ggplot(Expt123, aes(taxon, cMass)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic()

#pMass by taxon, split by experiment
ggplot(Expt123, aes(taxon, pMass)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic()

#pMass/hour by taxon, split by experiment
#first, general avg
Expt123$avg <- ((Expt123$pMass)/48) * 100

#then for expt3
sub3 <- which(Expt123$exptNo == 'Expt3')
Expt123[sub3,'avg'] <- (Expt123[sub3,'pMass']/168) *100
ggplot(Expt123, aes(taxon, avg)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic() + xlab('Taxon') + ylab('% Dissolution per Hour') +
  scale_fill_discrete(name = "Starting pH", labels = c("pH 5.1", "pH 7.1"))

#Work on this
stat_boxplot()
table <- A$stats


#just the two controls
#you could redo these looking at individual experiments rather than the average
subControl <- subset(Expt123, Expt123$taxon == 'Calcite'|Expt123$taxon == 'Tridacna')
ggplot(subControl, aes(taxon, avg)) +
  geom_boxplot(aes(fill = factor(exptNo))) +
  theme_classic() + xlab('Taxon') + ylab('% Dissolution per Hour') +
  scale_fill_discrete(name = "Starting pH", labels = c("pH 5.1", "pH 7.1"))


dev.off()

#3.3 comparison plots but standardized for pH
#subsets to add ph in to dataset, calculated separately.
#H+ is measured in molarity (mol/liter)
Expt1_3$deltaH <- 0.0000060473
Expt1_3[(Expt1_3$exptID == 'T1.2'),'deltaH'] <-0.0000042753
Expt1_3[(Expt1_3$exptID == 'T1.3'),'deltaH'] <-0.0000043809
Expt1_3[(Expt1_3$exptID == 'T4.1'),'deltaH'] <-0.0000000262
Expt1_3[(Expt1_3$exptID == 'T4.2'),'deltaH'] <-0.0000000488
Expt1_3[(Expt1_3$exptID == 'T4.3'),'deltaH'] <-0.0000000832

#change mol/liter hydrogen to mol CaCO3
Expt1_3$predictCarb <- ((Expt1_3$deltaH * 30.1)/2)
sub3 <- which(Expt1_3$exptNo == 'Expt3')
Expt1_3[sub3,'predictCarb'] <- (Expt1_3[sub3,'deltaH'] *30.075)/2

#total carbonate per replicate in mg
cTotal <- aggregate(Expt1_3$cMass, by = list(Expt1_3$exptID),FUN=sum)
colnames(cTotal)<-c('exptID','cTotal')
Expt123 <- merge(Expt1_3, cTotal, by = 'exptID')

#cMass is currently in mg, needs to be in mol like H+
#below equation goes from calcium carbonate solid mg to g to mol (divided by liters)
Expt1_3$cTotalMol <- (Expt1_3$cTotal)/(100.0869*1000)

#plot predicted carbonate over actual carbonate in mols
plot(predictCarb~cTotalMol, data= Expt1_3)
