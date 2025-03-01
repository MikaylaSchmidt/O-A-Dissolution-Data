###Hannah log script
setwd("C:/Users/micke/OneDrive/Desktop/Ch1 data/Hanna Logs Modified")
Expt1 <- read.csv('./Expt 1 Full.csv')
Expt1$duration <- (Expt1$Interval)*5
plot(pH.1.1 ~ duration, data=Expt1, ylim = c(5.0, 6.3), xlab= 'Duration (minutes)', ylab = 'Seawater pH', type='n')
log11 <- lm(pH.1.1 ~ duration, data=Expt1)
abline(log11 , col= 'pink')
log12 <- lm(pH.1.2 ~ duration, data=Expt1)
abline(log12, col='red')
log13 <- lm(pH.1.3 ~ duration, data=Expt1)
abline(log13, col='green')
legend('topleft', legend)


#now for expt 2
Expt2 <- read.csv('./Expt 2 Full.csv')
Expt2$duration <- (Expt2$Interval)*5
plot(pH.2.1 ~ duration, data=Expt2, ylim = c(5.0, 6.3), xlab= 'Duration (minutes)', ylab = 'Seawater pH', type='n')
log21 <- lm(pH.2.1 ~ duration, data=Expt2)
abline(log21 , col= 'pink')
log22 <- lm(pH.2.2 ~ duration, data=Expt2)
abline(log22, col='red')
log23 <- lm(pH.2.3 ~ duration, data=Expt2)
abline(log23, col='green')
legend('topleft', legend)

#now individual replicates
ExptValues <- data.frame(exptID = NA, meanTemp = NA, pH.01 = NA, pH.99 = NA, deltaPH = NA)

Expt1.1 <- read.csv('./Expt 1.1 Log 1.csv')
mean(Expt1.1$Temp...C.)
quantile(Expt1.1$pH, c(.01,.99))

Expt1.2 <- read.csv('./Expt 1.2 Log 1.csv')
mean(Expt1.2$Temp...C.)
quantile(Expt1.2$pH, c(.01,.99))

Expt1.3 <- read.csv('./Expt 1.3 Log 1.csv')
mean(Expt1.3$Temp...C.)
quantile(Expt1.3$pH, c(.01,.99))

Expt2.1 <- read.csv('./Expt 2.1 Log 1.csv')
mean(Expt2.1$Temp...C.)
quantile(Expt2.1$pH, c(.01,.99))

Expt2.2 <- read.csv('./Expt 2.2 Log 1.csv')
mean(Expt2.2$Temp...C.)
quantile(Expt2.2$pH, c(.01,.99))

Expt2.3 <- read.csv('./Expt 2.3 Log 1.csv')
mean(Expt2.3$Temp...C.)
quantile(Expt2.3$pH, c(.01,.99))

Expt3.1 <- read.csv('./Expt 3.1 Log Full.csv')
mean(Expt3.1$Temp...C.)
quantile(Expt3.1$pH, c(.01,.99))

Expt3.2 <- read.csv('./Expt 3.2 Log Full.csv')
mean(Expt3.2$Temp...C.)
quantile(Expt3.2$pH, c(.01,.99))


Expt3.3 <- read.csv('./Expt 3.3 Log Full.csv')
mean(Expt3.3$Temp...C.)
quantile(Expt3.3$pH, c(.01,.99))

