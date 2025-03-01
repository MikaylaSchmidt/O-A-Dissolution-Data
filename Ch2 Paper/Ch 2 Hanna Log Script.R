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
#Ch 2 Hanna Log Reading Script
#set working directory to ch2 hanna logs
setwd("C:/Users/micke/OneDrive/Desktop/Ch2 data/Ch 2 Hanna CSV")

Expt1.1 <- read.csv('./Ch2 Expt 1.1.csv')
mean(Expt1.1$Temp...C.)
quantile(Expt1.1$pH, c(.01,.99))

Expt1.2 <- read.csv('./Ch2 Expt 1.2.csv')
mean(Expt1.2$Temp...C.)
quantile(Expt1.2$pH, c(.01,.99))

Expt1.3 <- read.csv('./Ch2 Expt 1.3.csv')
mean(Expt1.3$Temp...C.)
quantile(Expt1.3$pH, c(.01,.99))

Expt2.1 <- read.csv('./Ch2 Expt 2.1.csv')
mean(Expt2.1$Temp...C.)
quantile(Expt2.1$pH, c(.01,.99))

Expt2.2 <- read.csv('./Ch2 Expt 2.2.csv')
mean(Expt2.2$Temp...C.)
quantile(Expt2.2$pH, c(.01,.99))

Expt2.3 <- read.csv('./Ch2 Expt 2.3.csv')
mean(Expt2.3$Temp...C.)
quantile(Expt2.3$pH, c(.01,.99))

Expt3.1 <- read.csv('./Ch2 Expt 3.1.csv')
mean(Expt3.1$Temp...C.)
quantile(Expt3.1$pH, c(.01,.99))

Expt3.2 <- read.csv('./Ch2 Expt 3.2.csv')
mean(Expt3.2$Temp...C.)
quantile(Expt3.2$pH, c(.01,.99))

Expt3.3 <- read.csv('./Ch2 Expt 3.3.csv')
mean(Expt3.3$Temp...C.)
quantile(Expt3.3$pH, c(.01,.99))

