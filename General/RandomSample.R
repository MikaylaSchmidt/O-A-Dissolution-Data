library(readxl)
Random_Sample <- read_excel("C:/Users/micke/OneDrive/Desktop/Random Sample.xlsx")
View(Random_Sample)
Block1 <- subset(Random_Sample, Block =="1")
View(Block1)
sample(Block1, size =12)
#ignore the block sample value, can't remember how to get rid of that
Block2 <- subset(Random_Sample, Block =="2")
View(Block2)
sample(Block2, size =12)
Block3 <- subset(Random_Sample, Block =="3")
View(Block3)
sample(Block3, size =12)
Block4 <- subset(Random_Sample, Block =="4")
View(Block4)
sample(Block4, size =12)
Block5 <- subset(Random_Sample, Block =="5")
View(Block5)
sample(Block5, size =12)

Random_Sample <- read.csv("C:/Users/micke/OneDrive/Desktop/Random Sample 2.1.csv")
View(Random_Sample)
sample(Random_Sample, size=9)
