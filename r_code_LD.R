pd<-read.csv(file.choose())
#subset the variables (subst) needed for analysis
vars_interest<-c("Group","Gender","Age","HYstage0", "UPDRS0","UPDRS4","UPDRS10","UPDRS16","LEDD0","LEDD4","LEDD10","LEDD16")
subst<-pd[vars_interest]

#look at summary to conduct a preliminary examination
prelim<-summary(subst)

#calculate averages of the LEDD and the UPRDS across 4 to 16 months
subst$UPDRSa <- rowMeans(subst[,6:8],na.rm = TRUE)
subst$LEDDa <- rowMeans(subst[,10:12],na.rm = TRUE)

#clean dataframe and select only those without NA's from dataframe subst in new data frame newpd
new_vars<-c("Group","Gender","HYstage0", "UPDRS0","UPDRSa","LEDD0","LEDDa")
temppd<-subst[new_vars]
newpd <- na.omit(temppd)

#calculate the changes in the outcome variable (LEDDa-LEDD0) and the mediator variable(UPDRSa-UPDRS0)
newpd$LEDDchange<-(newpd$LEDDa-newpd$LEDD0)
newpd$UPDRSchange<-(newpd$UPDRSa-newpd$UPDRS0)

#recode original group variable (FBF=5 AE=6 Ctrl=4) into newgr variable (FBF=2 AE=1 Ctrl=0)
newpd$newgr<-0
newpd[newpd$Group==5,]$newgr<-2
newpd[newpd$Group==6,]$newgr<-1

#test linear models; gender, hystage0 are control variables
#Effect of Predictor newgr on outcome 
fit1 <- lm(LEDDchange ~ Gender + newgr + HYstage0, data=newpd)
#Adjusted effect of the Predictor on Outcome
fit2 <- lm(LEDDchange ~ Gender + newgr + HYstage0+UPDRSchange, data=newpd)
#Effect of the predictor on the Mediator
fit3 <- lm(UPDRSchange ~ Gender + newgr + HYstage0, data=newpd)

#get standard error from covariance matrix for linear models
cov1<-sqrt(diag(vcov(fit1)))
cov2<-sqrt(diag(vcov(fit2)))
cov3<-sqrt(diag(vcov(fit3)))

#calculate indirect effect
indirectE<-fit1$coefficients[[3]]-fit2$coefficients[[3]]
#calculate standard error of direct effect using covariance matrices
gammax<-fit3$coefficients[[3]]
betaz<-fit2$coefficients[[5]]
segx<-cov3[[3]]
sebz<-cov2[[5]]
seIndirect<-sqrt((gammax^2)*(sebz^2)+(betaz^2)*(segx^2))
zscore<-indirectE/seIndirect
pval1<-pnorm(zscore)*2

#calculate confidence interval
my.ci.1<-indirectE-1.96*(seIndirect)
my.ci.u<-indirectE+1.96*(seIndirect)

#calculate proportion mediated
pmed<-indirectE*100/fit1$coefficients[[3]]
# confidence interval proportion mediated
pm.ci.l<-my.ci.1*100/fit1$coefficients[[3]] #actual upper limit, sign flip
pm.ci.u<-my.ci.u*100/fit1$coefficients[[3]] #actual lowerlimit

#predicted lms were not included in the final write-up, but they are here for reference

#plot lm , although not included in anaysis
# save predictions of the model in the new data frame 
# together with variable you want to plot against
library(ggplot2)
predicted_fit1 <- data.frame(pred1 = predict(fit1, newpd),newgr=newpd$newgr)

# this is the predicted line of multiple linear regression
ggplot(data = newpd, aes(x = newgr , y = LEDDchange)) + 
  geom_point(color='blue') +
  geom_line(color='red',data = predicted_fit1, aes(x=newgr, y=pred1))
#repeat above for adjusted regression
predicted_fit2 <- data.frame(pred2 = predict(fit2, newpd),newgr=newpd$newgr)

# this is the predicted line of multiple linear regression
ggplot(data = newpd, aes(x = newgr , y = LEDDchange)) + 
  geom_point(color='blue') +
  geom_line(color='red',data = predicted_fit2, aes(x=newgr, y=pred2))

library(plotly)

# Exploring bootstrap sampling distribution of LEDDchange across all 3 groups with bar plots
set.seed(515)
B <- 1000
my.boot <- numeric(B)
for (i in 1:B){
  x <- sample(newpd$LEDDchange, size=106, replace=TRUE)#draw resample that matches original sample size
  my.boot[i] <- mean(x) #compute mean, store in my.boot
}
hist(my.boot,xlab=expression(bar(X)),ylab="Density", main="Bootstrap distribution")
boot.mean <- mean(my.boot)
lines(c(boot.mean, boot.mean), c(0,250), col="blue", lwd=2, lty=2)

# check sample size for each group (0,1,2) for bootstrap resampling size
table(newpd$newgr)

#bootstrap for group 0 control
set.seed(515)
my.boot2 <- numeric(B)
for (i in 1:B){
  x2 <- sample(newpd$LEDDchange[newpd$newgr==0], size=36, replace=TRUE)#draw resample 
  my.boot2[i] <- mean(x2) #compute mean, store in my.boot
}
hist(my.boot2,xlab=expression(bar("Control LEDC")),ylab="Density", main="Bootstrap distribution")
boot.mean2 <- mean(my.boot2) #bootstrap mean for control
lines(c(boot.mean2, boot.mean2), c(0,250), col="blue", lwd=2, lty=2)

#Obtain Normal percentile 95% CI of control
controlLL <- boot.mean2-1.96*sd(my.boot2) #Lower limit of 95% Normal CI
controlUL <- boot.mean2+1.96*sd(my.boot2) #Upper limit of 95% Normal CI


#Regular ex group 1 or AE
set.seed(515)
my.boot3 <- numeric(B)
for (i in 1:B){
  x3 <- sample(newpd$LEDDchange[newpd$newgr==1], size=34, replace=TRUE)#draw resample 
  my.boot3[i] <- mean(x3) #compute mean, store in my.boot
}
hist(my.boot3,xlab=expression(bar("AE LEDC")),ylab="Density", main="Bootstrap distribution")
boot.mean3 <- mean(my.boot3)
lines(c(boot.mean3, boot.mean3), c(0,250), col="blue", lwd=2, lty=2)
#Obtain Normal percentile 95% CI of AE
aeLL <- boot.mean3-1.96*sd(my.boot3) #Lower limit of 95% Normal CI
aeUL <- boot.mean3+1.96*sd(my.boot3) #Upper limit of 95% Normal CI


#FBF group 2 bootstrap 
set.seed(515)
my.boot4 <- numeric(B)
for (i in 1:B){
  x4 <- sample(newpd$LEDDchange[newpd$newgr==2], size=36, replace=TRUE)#draw resample 
  my.boot4[i] <- mean(x4) #compute mean, store in my.boot
}
hist(my.boot4,xlab=expression(bar("FBF LEDC")),ylab="Density", main="Bootstrap distribution")
boot.mean4 <- mean(my.boot4)
lines(c(boot.mean4, boot.mean4), c(0,250), col="blue", lwd=2, lty=2)
#Obtain Normal percentile 95% CI of FBF
fbfLL <- boot.mean4-1.96*sd(my.boot4) #Lower limit of 95% Normal CI
fbfUL <- boot.mean4+1.96*sd(my.boot4) #Upper limit of 95% Normal CI


#summaries
summary(fit1)
summary(fit2)
summary(fit3)
