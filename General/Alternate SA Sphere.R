#alternate script for calculating ImageJ surface area of Eth and Nat

#calcSA1 = pi * r^2
dShell[subSphere,'r1'] <- sqrt(dShell[subSphere,'calcSA1'])/pi

#calcSA2 = pi * r^2
dShell[subSphere,'r2'] <- sqrt(dShell[subSphere,'calcSA2'])/pi

#SA is 4 * pi * r^2
#with the two different rs 
dShell[subSphere,'calcSA'] <- 4 * pi * (dShell[subSphere,'r1']*dShell[subSphere,'r2'])


#or, 
#without any of the calcSA1/2 values
dShell[subSphere,'calcSA'] <- 4 * pi * ((dShell[subSphere,'calcX']/2+dShell[subSphere,'calcY']/2 + dShell[subSphere,'yDim']/2)/3)^2
