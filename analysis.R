#################
# PREPARATIONS #
###############

# Load packages
library(plyr)
library(lme4)
library(pROC)
library(ggplot2)

# Read data into R
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_elisa/data.csv", header=TRUE, stringsAsFactors=FALSE)

# Format factor variables
data$Name <- as.factor(data$Name)
data$Age <- ordered(data$Age, levels=c("I", "SA", "A"))

# Format date variables
data$SmpDate <- as.POSIXct(data$SmpDate, format="%m/%d/%y")
data$DOB <- as.POSIXct(data$DOB, format="%m/%d/%y")

# Convert ODV values below threshold to 0, above to 1
data$Positive <- ifelse(data$ODV < 42.1, 0, 1)

########
# ROC #
######

# subset data
data_roc <- data[data$Age=="I" | data$Cyst==1,c("Cyst","ODV")]

# generate ROC curve data
ROC <- roc(data_roc$Cyst, data_roc$ODV)

# make data frame with roc curve data
df <- data.frame(TPR=ROC$sensitivities, FPR=1-ROC$specificities, threshold=ROC$thresholds)
df[df$threshold=="-Inf","threshold"] <- min(df$threshold)
df[df$threshold=="Inf","threshold"] <- max(df$threshold)

# plot ROC curve
qplot(df$FPR,df$TPR,geom="path",lwd=I(2),col=df$threshold,xlab="False positive rate",ylab="True positive rate")	+ scale_colour_gradient("Threshold") + geom_abline(slope=1,lty=2) + theme_classic()

##########################################################
# SUMMARY OF SAMPLES FROM INDIVIDUALS OF UNKNOWN STATUS #
########################################################

# subset data
data_unk <- data[data$Cyst==0 & data$Age != "I",]
data_unk$Age <- factor(data_unk$Age, ordered=F)

# summarize samples
nrow(data_unk) # 412 total rows of data
table(data_unk$Age) # 100 subadult, 312 adult sampes
table(data_unk$Male) # 202 female, 210 male samples

# summarize data at individual (not sample) level
indivs <- ddply(data_unk, .(Name), function(x) {
  Age <- x$Age[1]
  Male <- x$Male[1]
  Years <- mean(x$Years)
  numSamps <- nrow(x)
  numPos <- sum(x$Positive)
  allNeg <- ifelse(numPos==0, TRUE, FALSE)
  allPos <- ifelse(numSamps==numPos, TRUE, FALSE)
  trans <- ifelse(!allNeg & !allPos, TRUE, FALSE)
  data.frame(Age, Male, Years, numSamps, numPos, allNeg, allPos, trans)
  })

# demographic breakdown of individuals
nrow(indivs) # 158 unique individuals
table(indivs$Male) # 94 females, 64 males
table(indivs$Age) # 60 juvenile/subadult, 98 adult
range(indivs$Years) # age range 1.3-26.3 years

# distribution of positive samples
median(indivs$numSamps) # median number of samples per individual was 2
range(indivs$numSamps) # number of samples ranged from 1 to 14
sum(indivs$numPos) # 50 positive samples
sum(indivs$numSamps) - sum(indivs$numPos) # 362 negative samples
table(indivs$allNeg) # 26 individuals with at least 1 positive sample
indivs2 <- indivs[indivs$allNeg==FALSE,] # (subset to individuals with positive samples)
table(indivs2$Male) # 14 females and 12 males with positive samples
table(indivs2$Age) # 6 subadults, and 20 adults with positive samples
table(indivs2$trans) # 8 were always positive, 18 were 'transient'

# focus on non-cyst individuals sampled at least 5 times
indivs3 <- indivs[indivs$numSamps > 4,]
nrow(indivs3) # 23 well sampled individuals
nrow(indivs3[indivs3$allNeg==T,]) # 12 are all negative
nrow(indivs3[indivs3$numPos==1,]) # 7 had just one positive sample
indivs3[indivs3$numPos>1,] # inspect the remaining 4

# probability of observating >49 false positives in 412 samples
falsepos <- 1/57
sum(dbinom(50:412, size=412, prob=falsepos)) # p < 10^-25

##################################################################
# VISUALIZE OD VALUE DISTRIBUTIONS FOR CYST AND NON-CYST INDIVS #
################################################################

# Define bar heights and combine into matrix for barplot
breaks <- seq(0, 6, by=0.1)
cyst0 <- hist(log(data$ODV[data$Cyst==0] + 15), breaks=breaks, plot=F)
cyst1 <- hist(log(data$ODV[data$Cyst==1] + 15), breaks=breaks, plot=F)
cyst01 <- rbind(cyst0$counts, cyst1$counts)

# Stacked barplot
barplot(cyst01, space=0, legend.text=c("No Cyst", "Cyst"), col=c("dodgerblue3", "gray"), xlab="log(Index Value + 15)*10", ylab="Count (stacked)", ylim=c(0,60))
axis(side=1, at=(0:6)*10, labels=0:6)
abline(v=log(42.1 + 15)*10, lwd=2, lty=2)

###########
# MODELS #
#########

# Logistic regression predicting cysts
model.a1 <- glm(Cyst ~ Male * Years, data=indivs, family="binomial")
model.a2 <- glm(Cyst ~ Male + Years, data=indivs, family="binomial")
model.a3 <- glm(Cyst ~ Male, data=indivs, family="binomial")
model.a4 <- glm(Cyst ~ Years, data=indivs, family="binomial")
AIC(model.a1) # 72.427
AIC(model.a2) # 72.268
AIC(model.a3) # 83.782
AIC(model.a4) # 70.818
summary(model.a4) # positive effect of age (odds ratio = 1.21; p < 0.001)

# Remove infants and individuals with cysts for GLMMs predicting positive samples
data2 <- data[data$Cyst==0 & data$Age != "I",]

# Binomial GLMM with ordered categorial age, individual random effect
model.b1 <- glmer(Positive ~ Male * Age + (1|Name), data=data_unk, family="binomial") 
model.b2 <- glmer(Positive ~ Male + Age + (1|Name), data=data_unk, family="binomial")
model.b3 <- glmer(Positive ~ Male + (1|Name), data=data_unk, family="binomial")
model.b4 <- glmer(Positive ~ Age + (1|Name), data=data_unk, family="binomial")
model.b5 <- glmer(Positive ~ (1|Name), data=data_unk, family="binomial")
AIC(model.b1) # 225.8092
AIC(model.b2) # 224.0175
AIC(model.b3) # 222.0706
AIC(model.b4) # 222.051
AIC(model.b5) # 220.0918
summary(model.b5)

# Binomial GLMM with continuous age, individual random effect
model.c1 <- glmer(Positive ~ Male * Years + (1|Name), data=data_unk, family="binomial") 
model.c2 <- glmer(Positive ~ Male + Years + (1|Name), data=data_unk, family="binomial")
model.c3 <- glmer(Positive ~ Male + (1|Name), data=data_unk, family="binomial")
model.c4 <- glmer(Positive ~ Years + (1|Name), data=data_unk, family="binomial")
model.c5 <- glmer(Positive ~ (1|Name), data=data_unk, family="binomial")
AIC(model.c1) # 225.8593
AIC(model.c2) # 224.0287
AIC(model.c3) # 222.0706
AIC(model.c4) # 222.0439
AIC(model.c5) # 220.0918
summary(model.c5)

########
# END #
######




