## install pROC package
install.packages("pROC")

## load package
library(pROC); library(ggplot2)

## read data
data <- read.csv("~/Dropbox/Geladas/India's Code/ROCData.csv", stringsAsFactors=FALSE)
data$true_positive=NA
data[data$Age=="I","true_positive"]=0
data[data$Age=="A"&data$Cyst=="Y","true_positive"]=1

## generate ROC curve data
ROC=roc(data$true_positive,data$ODV)

## Area under the curve
ROC$auc

## make data frame with roc curve data
#true positive rate=sensitivity, false positive rate = (1- specificity), threshold = the threshold to call pos/neg
df=data.frame(TPR=ROC$sensitivities,FPR=1-ROC$specificities,threshold=ROC$thresholds)

## change -Inf and Inf threshold values to the min and max threshold numbers (for plotting):
df[df$threshold=="-Inf","threshold"]=min(df$threshold)
df[df$threshold=="Inf","threshold"]=max(df$threshold)

## plot the ROC curve
qplot(df$FPR,df$TPR,geom="path",lwd=I(2),col=df$threshold,xlab="false positive rate",ylab="true positive rate")+scale_colour_gradient("threshold")+geom_abline(slope=1,lty=2)+theme(legend.position=c(0.8,0.3))