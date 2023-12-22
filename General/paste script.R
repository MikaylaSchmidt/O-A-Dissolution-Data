##copy paste script
## getting rid of this
calcSA <- (pData$calcSA1 + pData$calcSA2) * (1/2) * 2


#SA for cube
dShell$calcSA <- calcSA1 + calcSA2 + 2*(dShell$calcY * dShell$zDim) + 2*(dShell$calcX * dShell$zDim)

#next lines for surface area of cylindrical specimens
# need to figure this out better for Liloa...
#equal to 2piR^2 + 2*h*pi*r
subCyl <- which(dShell$shape == 'cylinder')
dShell[subCyl,'calcSA'] <- calcSA1 + calcSA2 + dShell[subCyl,'zDim']*pi*(dShell[subCyl,'calcY'] + dShell[subCyl,'calcX'])/2

#next turbo, a hemisphere
#is pi*r^2 + 2*pi*r*h
subHemi <- which(dShell$shape == 'hemisphere')
dShell[subHemi,'calcSA'] <-  calcSA1 + 2 * pi * (dShell[subHemi,'calcX']/2 + dShell[subHemi,'calcY']/2)/2 * dShell[subHemi,'zDim'] 

#next lines for surface area of domed specimens
# is this correct?
subDome <- which(dShell$shape == '2dome')
dShell[subDome,'calcSA'] <- dShell[subDome,'zDim'] * calcSA1 + dShell[subDome,'zDim'] * calcSA2


#then, rough surface area for rhomboid Halimeda
#surface area is 2(xy/2) + 4(zc) where c = ((x/2)^2) + (y/2)^2)^1/2 
subRhom <- which(dShell$shape == 'rhombus')
dShell[subRhom,'calcC'] <- sqrt((dShell[subRhom,'calcX']/2)^2 + (dShell[subRhom,'calcY']/2)^2)
dShell[subRhom,'cSA1'] <- calcSA1 + calcSA2 + 4 * dShell[subRhom,'zDim'] * dShell[subRhom,'calcC']

#scaph and liloa need to be classified differently for this 
#for Monday... this is where you left off :)

#finally, nat and eth spherical calc
#is 4 * pi * [(x/2 + y/2 + z/2)/3]^2
subSphere <- which(dShell$shape == 'sphere')
dShell[subSphere,'cSA1'] <- 4 * pi * ((dShell[subSphere,'xDim']/2+dShell[subSphere,'yDim']/2 + dShell[subSphere,'yDim']/2)/3)^2

