#################
# PREPARATIONS #
###############

# Load packages
library(plyr)
library(lme4)

# Read data into R
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_elisa/data.csv", header=TRUE, stringsAsFactors=FALSE)

# Create factor variables
data$Name <- as.factor(data$Name)
data$Age <- ordered(data$Age, levels=c("I", "SA", "A"))

# Create date variables
data$SmpDate <- as.POSIXct(data$SmpDate, format="%m/%d/%y")
data$DOB <- as.POSIXct(data$DOB, format="%m/%d/%y")

# Convert ODV values below threshold to 0, above to 1
data$Positive <- ifelse(data$ODV < 42.1, 0, 1)

####################
# SUMMARY OF DATA #
##################

# demographic breakdown of samples
nrow(data) # 527 total rows of data
table(data$Male) # 261 females, 266 males
table(data$Age) # 57 infants, 104 juvenile/subadult, 366 adult

# summarize data for unique individuals
sum1 <- ddply(data, .(Name), function(x) {
  Age <- x$Age[1]
  Male <- x$Male[1]
  Years <- mean(x$Years)
  Cyst <- x$Cyst[1]
  numSamps <- nrow(x)
  numPos <- sum(x$Positive)
  allNeg <- ifelse(numPos==0, TRUE, FALSE)
  allPos <- ifelse(numSamps==numPos, TRUE, FALSE)
  trans <- ifelse(!allNeg & !allPos, TRUE, FALSE)
  data.frame(Age, Male, Years, Cyst, numSamps, numPos, allNeg, allPos, trans)
  })

# demographic breakdown of individuals
nrow(sum1) # 204 unique individuals
table(sum1$Male) # 117 females, 87 males
table(sum1$Age) # 37 infants, 60 juvenile/subadult, 107 adult
range(sum1$Years) # age range from ~0-26 years

# distribution of cysts 
sum(sum1$Cyst) # 10 cysts in 204 individuals (0.05 prevalence)
table(sum1[sum1$Cyst==1,"Male"]) # 6 females, 4 males
table(sum1[sum1$Cyst==1,"Age"]) # 0 infants, 1 subadult, 9 adults

# distribution of positive samples
median(sum1$numSamps) # median number of samples per individual was 2
range(sum1$numSamps) # number of samples ranged from 1 to 14
table(data$Positive) # 419 total negative samples, 108 total positive samples
table(sum1$allNeg) # 37 individuals with at least 1 positive sample
sum2 <- sum1[sum1$allNeg==FALSE,] # (subset to individuals with positive samples)
table(sum2$Male) # 20 females and 17 males with positive samples
table(sum2$Age) # 1 infant, 7 subadults, and 29 adults with positive samples
table(sum2$Cyst) # 10 cyst and 27 non-cyst individuals with positive samples
table(sum2$trans) # 17 were always positive, 20 were 'transient'
sum2[sum2$trans==TRUE,c("Name", "Cyst")] # only 1 cyst individual was 'transient' (Mary, 1 negative smp)

# focus on non-cyst individuals sampled at least 5 times
sum3 <- sum1[sum1$Cyst==0 & sum1$numSamps > 4,]
nrow(sum3) # 23 well sampled individuals
nrow(sum3[sum3$allNeg==T,]) # 12 are all negative
nrow(sum3[sum3$numPos==1,]) # 7 had just one positive sample
sum3[sum3$numPos>1,] # inspect the remaining 4

# probability of observating >50 false positives in 469 samples
falsepos <- 0.0179
1-sum(dbinom(0:50, size=469, prob=falsepos))

##################################################################
# VISUALIZE OD VALUE DISTRIBUTIONS FOR CYST AND NON-CYST INDIVS #
################################################################

# Define bar heights and combine into matrix for barplot
breaks <- seq(0, 6, by=0.1)
cyst0 <- hist(log(data$ODV[data$Cyst==0] + 15), breaks=breaks, plot=F)
cyst1 <- hist(log(data$ODV[data$Cyst==1] + 15), breaks=breaks, plot=F)
cyst01 <- rbind(cyst0$counts, cyst1$counts)

# Stacked barplot
quartz()
barplot(cyst01, space=0, legend.text=c("No Cyst", "Cyst"), col=c("gray28", "gray87"), xlab="log(Index Value + 15)*10", ylab="Count (stacked)", ylim=c(0,60))
axis(side=1, at=(0:6)*10, labels=0:6)
abline(v=log(42.1 + 15)*10, lwd=2, lty=2)

###########
# MODELS #
#########

# Logistic regression predicting cysts
model.a1 <- glm(Cyst ~ Male * Years, data=sum1, family="binomial")
model.a2 <- glm(Cyst ~ Male + Years, data=sum1, family="binomial")
model.a3 <- glm(Cyst ~ Male, data=sum1, family="binomial")
model.a4 <- glm(Cyst ~ Years, data=sum1, family="binomial")
AIC(model.a1) # 72.427
AIC(model.a2) # 72.268
AIC(model.a3) # 83.782
AIC(model.a4) # 70.818
summary(model.a4) # positive effect of age (odds ratio = 1.21; p < 0.001)

# Binomial GLMM with ordered categorial age, individual random effect
#model.b1 <- glmer(Positive ~ Male * Age + (1|Name), data=data, family="binomial") # Error: 'Downdated VtV is not positive definite'
model.b2 <- glmer(Positive ~ Male + Age + (1|Name), data=data, family="binomial")
model.b3 <- glmer(Positive ~ Male + (1|Name), data=data, family="binomial")
model.b4 <- glmer(Positive ~ Age + (1|Name), data=data, family="binomial")
model.b5 <- glmer(Positive ~ (1|Name), data=data, family="binomial")
#AIC(model.b1)
AIC(model.b2) # 283.7679
AIC(model.b3) # 280.6874
AIC(model.b4) # 281.7857
AIC(model.b5) # 278.6875
summary(model.b5)

# Binomial GLMM with continuous age, individual random effect
model.c1 <- glmer(Positive ~ Male * Years + (1|Name), data=data, family="binomial") 
model.c2 <- glmer(Positive ~ Male + Years + (1|Name), data=data, family="binomial")
model.c3 <- glmer(Positive ~ Male + (1|Name), data=data, family="binomial")
model.c4 <- glmer(Positive ~ Years + (1|Name), data=data, family="binomial")
model.c5 <- glmer(Positive ~ (1|Name), data=data, family="binomial")
AIC(model.c1) # 283.7513
AIC(model.c2) # 281.9724
AIC(model.c3) # 280.6874
AIC(model.c4) # 279.9943
AIC(model.c5) # 278.6875
summary(model.c5)

########
# END #
######




