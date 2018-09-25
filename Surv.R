library(car)
library(survival)
library(flexsurv)
library(KMsurv)
library(e1071)
library(rms)


gbcs <- read.csv("gbcs.csv", header = T)

names(gbcs)
## The data is from German Breast Cancer Study
## Research question is what influences the reoccurce of the breast cancer?
summary(gbcs)
hist(gbcs$age)
plot(density(gbcs$age))
table(gbcs$menopause) # 1 = no, 2 = yes
table(gbcs$hormone) # 1 = no, 2 = yes
hist(gbcs$prog_recp)
plot(density(gbcs$prog_recp), main = "Progesterone receptors, density plot")
# creating log transformation of progesterone receptor distribution
par(mfrow= c(1,2))
plot((gbcs$prog_recp), main = "Progesterone receptors, density plot")
plot((log(gbcs$prog_recp)), main = "Log of progesterone receptors, density plot")
## Following graphs will show some relationship between 
#  survival time and some explanatory variable
plot(gbcs$rectime~gbcs$age)
plot(gbcs$rectime~gbcs$menopause)
plot(gbcs$rectime~gbcs$hormone)
plot(gbcs$rectime~gbcs$prog_recp)
plot(gbcs$rectime~gbcs$estrg_recp)
plot(gbcs$rectime~gbcs$censrec) ## interesting

# correlation matrix
cor(gbcs)
scatterplotMatrix(gbcs)

# survival analysis begins with creating survival object
recsurv <- Surv(gbcs$rectime, gbcs$censrec)

# first thing we often do with survival data is to do KM curve
# technique to find out the proportion of people who are survivng over the time
fit_KM <- survfit(recsurv~1, type="kaplan-meier", conf.type = "log-log")
plot(fit_KM, main = "Survival function for rectime (K-M estimate)",
     xlab = "days", ylab="proportion")
# Other methods
fit_FH <- survfit(recsurv~1, type = "fleming-harrington", conf.type="log-log")
plot(fit_FH)

# print restricted means
print(fit_KM, print.rmean =T)

# survfits to illustrate impact of variables
leg.txt <- c("0", "1")
fit<-survfit(recsurv~as.numeric(gbcs$age>median(gbcs$age)))
plot(fit, col=c(2,4))
legend("topright", leg.txt, col=c(2,4), lty=1)
# with the above graph, we can observe that the 70% of younger patients 
# have survived without the reoccurence of the breast ca, but that numbers
# drops for older patients and therefore, older patients are more likely 
# to have the reoccurence of breast cancer

# for menopause
fit <- survfit(recsurv~gbcs$menopause)
plot(fit, col=c(1,3))
legend("topright", leg.txt,col=c(1,3), lty=1)

# for hormone therapy
fit <- survfit(recsurv~gbcs$hormone)
plot(fit, col=c(4,30), lty = c(1,2))
legend("topright", leg.txt,col=c(4,30), lty=c(1,2))

# for number of proges receptor
fit <- survfit(recsurv~as.numeric(gbcs$prog_recp>median(gbcs$prog_recp)))
plot(fit, col=c(1,2), lty = c(1,2))
legend("topright", leg.txt,col=c(1,2), lty=c(1,2))

# for number of estorgen receptor
fit <- survfit(recsurv~as.numeric(gbcs$estrg_recp>median(gbcs$estrg_recp)))
plot(fit, col=c(2,6), lty = c(1,1))
legend("topright", leg.txt,col=c(2,6), lty=c(1,1))

#