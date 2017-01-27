#################
# PREPARATIONS #
###############

# Load packages
library(plyr)
library(glmmADMB)

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
  checkCyst <- length(unique(x$Cyst))
  numSamps <- nrow(x)
  numPos <- sum(x$Positive)
  allNeg <- ifelse(numPos==0, TRUE, FALSE)
  allPos <- ifelse(numSamps==numPos, TRUE, FALSE)
  trans <- ifelse(!allNeg & !allPos, TRUE, FALSE)
  data.frame(Age, Male, Years, Cyst, checkCyst, numSamps, numPos, allNeg, allPos, trans)
  })

# Confirm that cyst status never changes across individuals during study period
unique(sum1$checkCyst)

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
sum2 <- sum1[sum1$allNeg==F,] # (subset to individuals with positive samples)
table(sum2$Male) # 20 females and 17 males with positive samples
table(sum2$Age) # 1 infant, 7 subadults, and 29 adults with positive samples
table(sum2$Cyst) # 10 cyst and 27 non-cyst individuals with positive samples
table(sum2$trans) # 17 were always positive, 20 were 'transient'
sum2[sum2$trans==TRUE,"Cyst"] # only 1 cyst individual was transient

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

# Which individual has a cyst and lies below the cutoff
data$Name[data$Cyst==1 & data$ODV<42.1] # Mary

###########
# MODELS #
#########

# Binomial GLM with ordered categorial age
model.b1 <- glm(Positive ~ Male + Age, data=data, family="binomial")
summary(model.b1) # positive linear effect of age (p < 0.01)

# Binomial GLM with continuous age
model.c1 <- glm(Positive ~ Male * Years, data=data, family="binomial")
summary(model.c1) # significant age and interaction term (age increases risk, stronger effect in males)

# Binomial GLMM with ordered categorial age, individual random effect
model.b2 <- glmmadmb(Positive ~ Male + Age + (1|Name), data=data, family="binom")
summary(model.b2) # nothing

# Binomial GLMM with continuous age, individual random effect
model.c2 <- glmmadmb(Positive ~ Male * Years + (1|Name), data=data, family="binom")
summary(model.c2) # nothing

########
# END #
######




