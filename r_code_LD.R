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

#plot lm's
plot(fit1)

# Exploring bootstrap sampling distribution of LEDDchange across all 3 groups
set.seed(515)
B <- 1000
my.boot <- numeric(B)
for (i in 1:B){
  x <- sample(newpd$LEDDchange, size=50, replace=TRUE)#draw resample 
  my.boot[i] <- mean(x) #compute mean, store in my.boot
}
hist(my.boot,xlab=expression(bar(X)),ylab="Density", main="Bootstrap distribution")
boot.mean <- mean(my.boot)
lines(c(boot.mean, boot.mean), c(0,250), col="blue", lwd=2, lty=2)
