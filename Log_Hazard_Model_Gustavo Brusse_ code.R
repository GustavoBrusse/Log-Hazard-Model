#### Gustavo Brusse - EDSD 2018/2019 ###
#### Assignment Event History 2 ########

## Exercise 1
## 1.1
library(survival) 

## Reading the dataset
setwd("C:\\Users\\gbrusse\\Desktop\\Event History 2")
nursing<- read.table("NursingHomeMod.txt",header=T)

## Creating a data frame
data<-cbind(ID=1:1601, nursing)

## split time axis and replicating the individuals
data2 <- survSplit(data, cut=c(0.25, 0.5, 1, 2, 3), end="stay", event="event", start="start", episode="interval")
data2$y.new<-data2$stay-data2$start
data2$interval<-as.factor(data2$interval)
head(data2)

## 1.2
## Model for the hazard without any further covariates
model1<-glm(event ~ interval, family = poisson, data = data2, offset = log(y.new)) 
summary(model1)

model2<-glm(event ~ -1, family = poisson, data = data2, offset = log(y.new)) 
summary(model2)

## 1.3 
## Ploting the results
taujj<-c(0,0.25,0.5,1,2,3)
haz1 <- model1$coeff[1]
haz2 <- model1$coeff[2]
haz3 <- model1$coeff[3]
haz4 <- model1$coeff[4]
haz5 <- model1$coeff[5]
haz<-c(haz1,haz2,haz3,haz4,haz5) # log Hazard
ht.hat <- exp(model1$coeff)
# Plot P-W Constant Hazard
plot(1, 1, t="n", ylim=c(0,3), xlim=c(0,3),
     ylab="h(t)", xlab="years",main="Piece-Wise Constant Hazard")
for(i in 1:length(taujj)){
  segments(x0=taujj[i], y0=ht.hat[i], x1=taujj[i+1], y1=ht.hat[i], lwd=6, col=1)
}

## 1.4
## Creating the models
ms<-glm(event ~ interval+sex, family = poisson, data = data2, offset = log(y.new)) 
summary(ms)
mSM<-glm(event ~ interval+sex+marstat, family = poisson, data = data2, offset = log(y.new))
summary(mSM)
mSMH<-glm(event ~ interval+sex+marstat+health, family = poisson, data = data2, offset = log(y.new))
summary(mSMH)
mSMHA<-glm(event ~ interval+sex+marstat+health+agegr, family = poisson, data = data2, offset = log(y.new))
summary(mSMHA)

## Comparing the models
anova(ms,mSM,mSMH,mSMHA, test="Chisq")

## Exercise 2
## reading the dataset
setwd("C:\\Users\\gbrusse\\Desktop\\Event History 2")
breast<-read.table("breast.txt",header=T)

## 2.1 
dati <- data.frame(id=1:686, exit=breast$exit, status=breast$event, 
                   group=breast$group)

## creating a X matrix
dati$status<-as.numeric(dati$status) 
good=ifelse(dati$group=="Good",1,0)
medium=ifelse(dati$group=="Medium",1,0)
poor=ifelse(dati$group=="Poor",1,0)
group<-matrix(c(good,medium),ncol = 2)

## likelyhood function
llkPHweibull <- function(par, exit, status, group){ 
  ## weibull parameters 
  a <- par[1] 
  b <- par[2] 
  ## regression coefficients 
  beta <- par[3:4] 
  ## design matrix 
  X <- group
  ## regression part 
  Xbeta <- X%*%beta 
  expXbeta <- exp(Xbeta) 
  ## log hazard - weibull a and b
  h0 <- log(a/b) + (a-1)*log(exit/b) 
  ## cumulative hazards 
  H0 <--log(1-pweibull(exit,a,b))
  ## individual contribution 
  llki <- status*(h0+Xbeta) - H0*expXbeta 
  ## total likelihood 
  LLK <-sum(llki) 
  return(LLK)
} 

## MLE Optimization
mle<-optim(c(1,10,1,1), llkPHweibull, exit=dati$exit,
           status=dati$status,group=group, control=list(fnscale=-1),hessian=T)
a.hat <- mle$par[1]
b.hat <- mle$par[2]
beta1.hat <- mle$par[3]
beta2.hat <- mle$par[4]

## 2.2
mle$hessian
vcov <- solve(-mle$hessian)
se.a <- sqrt(vcov[1,1])
se.b <- sqrt(vcov[2,2])
se.beta1 <- sqrt(vcov[3,3])
se.beta2 <- sqrt(vcov[4,4])

CIs <- data.frame(hat=c(a.hat,
                        b.hat,
                        beta1.hat,
                        beta2.hat),
                  low=c(a.hat-1.96*se.a,
                        b.hat-1.96*se.b,
                        beta1.hat-1.96*se.beta1,
                        beta2.hat-1.96*se.beta2),
                  up= c(a.hat+1.96*se.a,
                        b.hat+1.96*se.b,
                        beta1.hat+1.96*se.beta1,
                        beta2.hat+1.96*se.beta2))
rownames(CIs) <- c("a", "b", "beta1", "beta2")
CIs

# 2.3 Ploting the log-hazard 
exit <- 0:100
lh0 <-log(a.hat/b.hat)+log(exit/b.hat)*(a.hat-1)+ beta1.hat*0+beta2.hat*0
lh1 <-log(a.hat/b.hat)+log(exit/b.hat)*(a.hat-1)+ beta1.hat*0+beta2.hat*1
lh2 <-log(a.hat/b.hat)+log(exit/b.hat)*(a.hat-1)+ beta1.hat*1+beta2.hat*0
plot(exit, lh0, type="l",col=2,
     xlab="time", ylab="log-hazard",ylim=c(-5,2), main="Log-hazard")
lines(exit, lh1,col=4)
lines(exit, lh2,col=6)
legend("topleft", legend=c("Poor", "Medium","Good"),
       col=c(2,4,6),lty=c(1,1))

## Exercise 3
library(survival)
library(MASS)
china<- read.table("ChinaMod.txt",header=T)

china$gender<-as.factor(china$gender)
china$urbrur<-as.factor(china$urbrur)
china$act<-as.factor(china$act)
china$adl<-as.factor(china$adl)

##3.1
surv.dat <- Surv(time=china$entry, time2=china$exit, event=china$status) 
m1 <- coxph(surv.dat ~ gender*urbrur, data=china)
summary(m1)

## 3.2
full.cox <- coxph(surv.dat ~ gender*urbrur*act*adl, data=china)
step.cox<-stepAIC(full.cox) 
summary(step.cox) 

## 3.3
## Best model
best <- coxph(surv.dat ~ gender+urbrur+act+adl+urbrur*act+gender*adl, data=china)
summary(best)

## male, living in rural area, no limitations
male<- exp(best$coefficients[1]+best$coefficients[2])

## male, living in rural area, no limitations
female<- exp(best$coefficients[2])

## 
female/male

 
