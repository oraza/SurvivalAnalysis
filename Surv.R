library(car)
library(survival)
library(flexsurv)
library(KMsurv)
library(e1071)
library(rms)


gbcs <- read.csv("gbcs.csv", header = T)
# or get it from online 
# gbcs <- read.csv("https://ryanwomack.com/data/gbcs.csv")
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

##### survreg to check the distribution
## here we will identify the best fit
## with survreg function we will run regression with various families
## (distributions) on a survival object
fit.exp <- survreg(recsurv~1, dist = "exponential")
fit.wbull <- survreg(recsurv~1, dist = "weibull")
fit.gaus <- survreg(recsurv~1, dist = "gaussian")
fit.logis <- survreg(recsurv~1, dist = "logistic")
fit.lognrom <- survreg(recsurv~1, dist = "lognormal")
fit.loglogis <- survreg(recsurv~1, dist = "loglogistic")
summary(fit.exp)
summary(fit.wbull)
summary(fit.gaus)
summary(fit.logis)
summary(fit.lognrom)
summary(fit.loglogis)

## the same can be done with flexsurvreg function 
fit.exp.flex <- flexsurvreg(recsurv~1, dist = "exp")
fit.wbull.flex <- flexsurvreg(recsurv~1, dist = "weibull")
fit.gamma.flex <- flexsurvreg(recsurv~1, dist = "gamma")
fit.ggamma.flex <- flexsurvreg(recsurv~1, dist = "gengamma")
fit.genf.flex <- flexsurvreg(recsurv~1, dist = "genf")
fit.lognorm.flex <- flexsurvreg(recsurv~1, dist = "lnorm")
fit.gompertz.flex <- flexsurvreg(recsurv~1, dist = "gompertz")
fit.exp.flex 
fit.wbull.flex
fit.gamma.flex 
fit.ggamma.flex 
fit.genf.flex 
fit.lognorm.flex 
fit.gompertz.flex

# Plotting the fit
plot(fit.exp.flex)
plot(fit.wbull.flex)
plot(fit.gamma.flex)
plot(fit.ggamma.flex) 
plot(fit.genf.flex) 
plot(fit.lognorm.flex) 
plot(fit.gompertz.flex)

# extracting log likelood test values
fit.exp.flex$loglik
fit.wbull.flex$loglik
fit.gamma.flex$loglik
fit.ggamma.flex$loglik
fit.genf.flex$loglik
fit.lognorm.flex$loglik
fit.gompertz.flex$loglik

# significance test
sign_lognorm <- 2*(fit.ggamma.flex$loglik - fit.lognorm.flex$loglik)
sign_lognorm ## not significantly different from each other

## reporting AIC
fit.exp.flex$AIC
fit.wbull.flex$AIC
fit.gamma.flex$AIC
fit.ggamma.flex$AIC
fit.genf.flex$AIC
fit.lognorm.flex$AIC
fit.gompertz.flex$AIC

## re-running analysis with log transform 
logprog <- log(gbcs$prog_recp +0.1)
logest <- log(gbcs$estrg_recp +0.1)
logrectime <- log(gbcs$rectime)
logrecsurv <- Surv(logrectime, gbcs$censrec)

# finding a best fit with logrecsurv
fit.exp.log <- flexsurvreg(logrecsurv~1, dist = "exp")
fit.wbull.log <- flexsurvreg(logrecsurv~1, dist = "weibull")
fit.gamma.log <- flexsurvreg(logrecsurv~1, dist = "gamma")
fit.ggamma.log <- flexsurvreg(logrecsurv~1, dist = "gengamma")
fit.genf.log <- flexsurvreg(logrecsurv~1, dist = "genf")
fit.lognorm.log <- flexsurvreg(logrecsurv~1, dist = "lognormal")
fit.exp.log
fit.wbull.log
fit.gamma.log
fit.ggamma.log
fit.genf.log 
fit.lognorm.log 
# plotting
plot(fit.exp.log)
plot(fit.wbull.log)
plot(fit.gamma.log)
plot(fit.ggamma.log)
plot(fit.genf.log) 
plot(fit.lognorm.log) 

# loglikelihood test
fit.exp.log$loglik
fit.wbull.log$loglik
fit.gamma.log$loglik
fit.ggamma.log$loglik
fit.genf.log$loglik
fit.lognorm.log$loglik 

# AIC based
fit.exp.log$AIC
fit.wbull.log$AIC
fit.gamma.log$AIC
fit.ggamma.log$AIC
fit.genf.log$AIC
fit.lognorm.log$AIC 

# creating graphs
par(mfrow = c(1,2))
plot(fit.ggamma.flex, main = "Rectime fit to generalized gamma",
     xlab= "days", ylab = "p")
plot(fit.lognorm.log, main = " (Log) rectime fit to lognormal",
     xlab = "days", ylab = "p")

# CI and diffplots
# diffplot will show a differential survival within two groups
# to evaluate the observed difference by comparing the confidence
# interval around the observed survival curves
# we will use 'npsurv' funtion from rms library 
fit <- npsurv(recsurv~gbcs$hormone)
survplot(fit)

# quantitative way to assess this is by:
survdiffplot(fit)

# numerical test of whether variable influences survival
# rho = 0 logrank (Mantel-Haenzel) test
# rho = 1 Peto & peto modification of the Gehan-Wilcoxon test

survdiff(recsurv~ gbcs$hormone, rho = 0) 
survdiff(recsurv~ gbcs$menopause, rho = 0)
gbcs$agecat <- as.numeric(gbcs$age>median(gbcs$age))
survdiff(recsurv~ gbcs$agecat, rho = 0)
gbcs$prog_recp.cat <- as.numeric(gbcs$prog_recp>median(gbcs$prog_recp))
survdiff(recsurv~ gbcs$prog_recp.cat, rho = 0)
gbcs$estrg_recp.cat <- as.numeric(gbcs$estrg_recp>median(gbcs$estrg_recp))
survdiff(recsurv~ gbcs$estrg_recp.cat, rho = 0)
# Repeating with Peto & Peto modification
survdiff(recsurv~ gbcs$hormone, rho = 1) #
survdiff(recsurv~ gbcs$menopause, rho = 1)
survdiff(recsurv~ gbcs$agecat, rho = 1)
survdiff(recsurv~ gbcs$prog_recp.cat, rho = 1)
survdiff(recsurv~ gbcs$estrg_recp.cat, rho = 1)

# Cox proportional hazard model with coxph and testing fit
fit.cox1 <- coxph(recsurv~gbcs$age)
fit.cox1
cox.zph(fit.cox1)
plot(cox.zph(fit.cox1))

fit.cox2 <- coxph(recsurv~gbcs$menopause)
fit.cox2
cox.zph(fit.cox2)
plot(cox.zph(fit.cox2))

fit.cox3 <- coxph(recsurv~gbcs$hormone)
fit.cox3
cox.zph(fit.cox3)
plot(cox.zph(fit.cox3))

fit.cox4 <- coxph(recsurv~gbcs$prog_recp)
fit.cox4
cox.zph(fit.cox4)
plot(cox.zph(fit.cox4))

fit.cox5 <- coxph(recsurv~gbcs$estrg_recp)
fit.cox5
cox.zph(fit.cox5)
plot(cox.zph(fit.cox5))