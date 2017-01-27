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

#########################################################################
# VISUALIZE POSITIVE / NEGATIVE SAMPLES
#########################################################################



#########################################################################
# VISUALIZE OD VALUE DISTRIBUTIONS FOR CYST AND NON-CYST INDIVS
#########################################################################

# Define bar heights and combine into matrix for barplot
breaks <- seq(0, 6, by=0.1)
cyst0 <- hist(log(data$ODV[data$Cyst==0] + 14), breaks=breaks, plot=F)
cyst1 <- hist(log(data$ODV[data$Cyst==1] + 14), breaks=breaks, plot=F)
cyst01 <- rbind(cyst0$counts, cyst1$counts)

# Stacked barplot
quartz()
barplot(cyst01, space=0, legend.text=c("No Cyst", "Cyst"), col=c("gray28", "gray87"), xlab="log(Index Value + 14)", ylab="Count (stacked)", ylim=c(0,60))
axis(side=1, at=(0:6)*10, labels=0:6)
abline(v=log(42.1 + 14)*10, lwd=2, lty=2)

# Which individual has a cyst and lies below the cutoff
data$Name[data$Cyst==1 & data$ODV<37.5] # Mary

#########################################################################
# VISUALIZE SAMPLE Positive AND DATES FOR ALL SK INDIVIDUALS
#########################################################################

# Create data for figures
fig1 <- data[data$Name %in% x2$Name,]
fig1m <- fig1[fig1$Sex == "M",]
fig1f <- fig1[fig1$Sex == "F",]

# Names ordered by age
m.names <- ddply(fig1m, .(Name), function(x){
  x <- x[order(x$SmpDate),]
  data.frame(x$DOB[1], x$Age[1])
})
f.names <- ddply(fig1f, .(Name), function(x){
  x <- x[order(x$SmpDate),]
  data.frame(x$DOB[1], x$Age[1])
})
m.names <- m.names[order(m.names[,2], decreasing=T),]
m.names[,3] <- as.character(m.names[,3])
f.names <- f.names[order(f.names[,2], decreasing=T),]
f.names[,3] <- as.character(f.names[,3])
m.names[m.names[,3]=="I",3] <- paste(m.names[m.names[,3]=="I",3], 1:sum(m.names[,3]=="I"), sep="")
m.names[m.names[,3]=="J",3] <- paste(m.names[m.names[,3]=="J",3], 1:sum(m.names[,3]=="J"), sep="")
m.names[m.names[,3]=="SA",3] <- paste(m.names[m.names[,3]=="SA",3], 1:sum(m.names[,3]=="SA"), sep="")
m.names[m.names[,3]=="A",3] <- paste(m.names[m.names[,3]=="A",3], 1:sum(m.names[,3]=="A"), sep="")
f.names[f.names[,3]=="I",3] <- paste(f.names[f.names[,3]=="I",3], 1:sum(f.names[,3]=="I"), sep="")
f.names[f.names[,3]=="J",3] <- paste(f.names[f.names[,3]=="J",3], 1:sum(f.names[,3]=="J"), sep="")
f.names[f.names[,3]=="SA",3] <- paste(f.names[f.names[,3]=="SA",3], 1:sum(f.names[,3]=="SA"), sep="")
f.names[f.names[,3]=="A",3] <- paste(f.names[f.names[,3]=="A",3], 1:sum(f.names[,3]=="A"), sep="")

# Two panel plot
quartz()
layout(matrix(1:2, 1, 2))

# Males
plot(c(1:10)~c(1:10), ylim=c(0,nrow(m.names)), xlim=c(min(fig1m$SmpDate), max(fig1m$SmpDate)), type="n", xaxt="n", yaxt="n", xlab="Sample date", ylab="Individual", main="Males (n = 50)")
ddply(fig1[fig1$Sex=="M",], .(Name), function(x) {
  x <- x[!duplicated(x$SmpDate),]
  x <- x[order(x$SmpDate),]
  cols <- ifelse(x$Positive==1, "black", "white")
  pos <- which(m.names[,1]==x$Name[1])
  #points(rep(pos, nrow(x)) ~ SmpDate, data=x, type="l", col="gray")
  abline(h=pos, col="gray", lty=3)
  points(rep(pos, nrow(x)) ~ SmpDate, data=x, pch=21, bg=cols, cex=0.5)
})
axis(side=2, at=1:nrow(m.names), labels=m.names[,3], las=2, cex.axis=0.3, tcl=-0.2, mgp=c(3,0.3,0))
xdates <- c("Aug 2014", "Sep 2014", "Oct 2014", "Apr 2015", "May 2015", "June 2015")
xpos <- as.POSIXct(c("2014/08/01", "2014/09/01", "2014/10/01", "2015/04/01", "2015/05/01", "2015/06/01"), format="%Y/%m/%d")
axis(side=1, at=xpos, labels=FALSE)
text(x=xpos, y=rep(-7, length(xpos)), labels=xdates, cex=0.75, srt=45, pos=1, xpd=T)

# Add legend
legend(1402047062, 63.27607, legend=c("Positive", "Negative"), pch=21, pt.bg=c("black", "white"), xpd=TRUE, cex=0.75)

# Females
plot(c(1:10)~c(1:10), ylim=c(0,nrow(f.names)), xlim=c(min(fig1f$SmpDate), max(fig1f$SmpDate)), type="n", xaxt="n", yaxt="n", xlab="Sample date", ylab="", main="Females (n = 57)")
ddply(fig1[fig1$Sex=="F",], .(Name), function(x) {
  x <- x[!duplicated(x$SmpDate),]
  x <- x[order(x$SmpDate),]
  cols <- ifelse(x$Positive==1, "black", "white")
  pos <- which(f.names[,1]==x$Name[1])
  #points(rep(pos, nrow(x)) ~ SmpDate, data=x, type="l", col="gray")
  abline(h=pos, col="gray", lty=3)
  points(rep(pos, nrow(x)) ~ SmpDate, data=x, pch=21, bg=cols, cex=0.5)
})
axis(side=2, at=1:nrow(f.names), labels=f.names[,3], las=2, cex.axis=0.3, tcl=-0.2, mgp=c(3,0.3,0))
axis(side=1, at=xpos, labels=FALSE)
text(x=xpos, y=rep(-7, length(xpos)), labels=xdates, cex=0.75, srt=45, pos=1, xpd=T)

#########################################################################
# VISUALIZE INFECTION RATE BY AGE AND SEX
#########################################################################

# Bar heights for infection
indivs <- table(x3$Sex, x3$Age)[,c(2,3,4,1)]
indivs <- cbind(indivs, rowSums(indivs))
positive <- table(x3[x3$Pos==1,"Sex"], x3[x3$Pos==1,"Age"])[,c(2,3,4,1)]
positive <- cbind(positive, rowSums(positive))
prop.pos <- positive/indivs
colnames(prop.pos) <- colnames(positive) <- colnames(indivs) <- c("Infant", "Juvenile", "Subadult", "Adult", "Total")
err.upper <- c()
err.lower <- c()
for (i in 1:10) {
  test <- prop.test(positive[i], indivs[i])
  err.upper[i] <- test$conf.int[2]
  err.lower[i] <- test$conf.int[1]
}

# Barplot with error bars for proportion with at least one positive
quartz()
pp.bplot <- barplot(prop.pos, ylim=c(0, 1), legend.text=c("Female", "Male"), axis.lty=1, beside=T, ylab="Proportion infected", xlab="Age category")
arrows(x0=c(pp.bplot), y0=err.lower, x1=c(pp.bplot), y1=err.upper, length=0.07, angle=90, code=3)

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




