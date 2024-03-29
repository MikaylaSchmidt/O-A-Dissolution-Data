#means file
pData <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)
pData$deviation <- (abs(pData$xDim - pData$cSize) + abs(pData$yDim - pData$cSize) +  abs(pData$zDim - pData$cSize))/3
Expt123 <- pData

abranda <- subset(Expt123, Expt123$taxon =='Abranda')
summary(abranda)

calcite <- subset(Expt123, Expt123$taxon =='Calcite')
summary(calcite)

eth <- subset(Expt123, Expt123$taxon =='Ethalia')
summary(eth)

hali <- subset(Expt123, Expt123$taxon =='Halimeda')
summary(hali)

liloa <- subset(Expt123, Expt123$taxon =='Liloa')
summary(liloa)

margi <- subset(Expt123, Expt123$taxon =='Marginopora')
summary(margi)

natica <- subset(Expt123, Expt123$taxon =='Natica')
summary(natica)

pingui <- subset(Expt123, Expt123$taxon =='Pinguitellina')
summary(pingui)

scaph <- subset(Expt123, Expt123$taxon =='Scaphopod')
summary(scaph)

tridacna <- subset(Expt123, Expt123$taxon =='Tridacna')
summary(tridacna)

turbo <- subset(Expt123, Expt123$taxon =='Turbo')
summary(turbo)

