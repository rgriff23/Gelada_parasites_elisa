#################
# PREPARATIONS #
###############

# Load packages
library(plyr) # data wrangling
library(ggplot2) # plots
library(pROC) # ROC analysis
library(lme4) # GLMMs
library(MuMIn) # model comparison


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

# Subset to subadults/adults for regression predicting cysts
data_cysts <- ddply(data, .(Name), function(x) {data.frame(Cyst=x$Cyst[1], Male=x$Male[1], Years=mean(x$Years))})

# Logistic regression predicting cysts
mod1 <- glm(Cyst ~ Male * Years, data=data_cysts, family="binomial", na.action=na.fail)
dredge(mod1)

# Binomial GLMM with ordered categorial age, individual random effect
mod2 <- glmer(Positive ~ Male * Age + (1|Name), data=data_unk, family="binomial", na.action=na.fail) 
dredge(mod2)

# Binomial GLMM with continuous age, individual random effect
mod3 <- glmer(Positive ~ Male * Years + (1|Name), data=data_unk, family="binomial", na.action=na.fail) 
dredge(mod3)

########
# END #
######




