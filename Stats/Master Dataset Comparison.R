################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2024 Mikayla Schmidt
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

#Script for creating graphics to visualize combined data from dissolution Master Dataset
#Portions of code adapted from 'Expt2 Post Processing Analysis'

#0.Open clean master dataset 
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data")
dShell <- read.csv('./shellsData_Ch1MasterSet.csv', skip=26)

pData <- dShell
TAXA <- unique(pData$taxon)
pData$taxon <- as.factor(pData$taxon)
taxa <- sort(unique(pData$taxon))
tFont <- rep(3,length(taxa))
tFont[which(taxa == 'Abranda')] <- 1
taxaAbrev <- substring(taxa,0,5)

#1.1 Combining Expt 1 and Expt 2 without Expt 4
#first, you need to subset it
  
  Expt1_2 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T2.1'|pData$exptID == 'T2.2'|pData$exptID == 'T2.3')

#1.2 Setup and execution of split plots
  pdf('./outFigs/1&2SplitWax.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  #1. firstly, some setup to subset
  library(ggplot2)
  factor(Expt1_2$waxYN)
  Expt1_2$waxYN <- as.factor(Expt1_2$waxYN)
  str(Expt1_2)
  
  #2. cMass by taxon, split by wax presence
  ggplot(Expt1_2, aes(taxon, cMass)) +
    geom_boxplot(aes(fill = factor(waxYN))) +
    theme_classic()
  
  #3 pMass by taxon, split by wax presence
  ggplot(Expt1_2, aes(taxon, pMass)) +
    geom_boxplot(aes(fill = factor(waxYN))) +
    theme_classic()
  
  dev.off()  
  
  
#2.1 Combining Expt 1 and Expt 4 without Expt 2, and then subsetting by experiment
  
  Expt1_4 <- subset(pData, pData$exptID == 'T1.1'|pData$exptID == 'T1.2'|pData$exptID == 'T1.3'|pData$exptID == 'T4.1'|pData$exptID == 'T4.2'|pData$exptID == 'T4.3')
  Expt1_4$exptNo <- 'Expt1'
  Expt1_4[(Expt1_4$exptID == 'T4.1'),'exptNo'] <-'Expt4'
  Expt1_4[(Expt1_4$exptID == 'T4.2'),'exptNo'] <-'Expt4'
  Expt1_4[(Expt1_4$exptID == 'T4.3'),'exptNo'] <-'Expt4'
  
#2.2 Setup and execution of comparison plots
  #only thing necessary for this is to split it into two more subsets by overall experiment 
  
  pdf('./outFigs/1&4SplitExpt.pdf', page='A4', height = 6 , width = 8)
  par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(8,4,1,0))
  
  
  #1. firstly, some setup to subset
  library(ggplot2)
  factor(Expt1_4$exptNo)
  Expt1_4$exptNo <- as.factor(Exp1_4$exptNo)
  str(Expt1_4)
  
  #2. cMass by taxon, split by wax presence
  ggplot(Expt1_4, aes(taxon, cMass)) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic()
  
  #3 pMass by taxon, split by wax presence
  ggplot(Expt1_4, aes(taxon, pMass)) +
    geom_boxplot(aes(fill = factor(exptNo))) +
    theme_classic()
  
  dev.off()  
  