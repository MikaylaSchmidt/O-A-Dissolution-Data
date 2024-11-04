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
