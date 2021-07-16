
### regression approach with norm, simulates regression model, not common distribution y, x_1,..x_p

### samples defined by location of predictors and  covariance (pairwise correlated) 
### regression defined by deviation from true parameters

### version includes a loop to modify - for example sample sizes!!!


library(MASS)
library(car)
library(pracma)




n <- 40 ## probsample
n1 <- 200  ## nonprob from target
n2 <- 400  ## nonprob from poluting

numcont <- 4 ### number continuous var - in the following !!!! means have to be dapted (twice)
numpred <- numcont

###propsample
########################

##covariance  predictors 
var0  <- rep(1,numcont)  ## variance predictors
cov0 <- 0.3  ## pairwise covariance

covm <- diag(var0-cov0,numcont)+ matrix(cov0,numcont,numcont)  ## second var= dependent variable

### location 
mu <- rep(1,numcont)

### true regression
trueint <- 1   ## intercept
truepar <- seq(1,numcont,1)  ## parameters covariates
truevect <- c(trueint,truepar)  ## true vector

respstdev <- 1  ## stddeviation response


################

## pollution sample
##############################
var2  <- rep(1,numcont)  ## variance predictors
cov2 <- 0.3  ## pairwise covariance

covm2 <- diag(var0-cov0,numcont)+ matrix(cov0,numcont,numcont)  ## second var= dependent variable

### location 
mu2 <- mu + c(1,0,2,1  ) #!!!! length numcont

###  regression in polluted sample
pollint <- 2 ##intercept

pollpar <- truepar+ c(-2,0,1,2 ) ## !!!! length numcont parameters covariates


respstdevpol <-2  ## stddeviation

#######################


 

#### loop probsamplesize starting from 30
numsimpss <- 5 
resultspss <- matrix(0,numsimpss,2)
numobs<- matrix(0,numsimpss,1)

for (lp in 1:numsimpss){
n <- 30+20*(lp-1)
numobs[lp,1]<-n

#########################################
#### loop number of simulated data sets
numsim <- 10
 
resultprob <- matrix (0,numsim,numpred+1)
resultext <- matrix (0,numsim,numpred+1)

rates <- matrix (0,numsim,2)

#repet <- 1 #number iterations within obsolete
#change <- matrix(0,numsim,repet )

threshold <- 1.96 ## for residuals
#critbeta <- 0.1   ## for relative change of beta, now computed, replaced by quantbetachange
quantbetachange <- 0.95 #quantile for beta change

for (l in 1: numsim) {
  
  ### prob sample  
  
  prob <- mvrnorm(n , mu, covm, tol = 1e-6, empirical = FALSE )
  #colnames(prob) <- c("X1", "X2", "X3")
  prob <- as.data.frame(prob)
  prob$lab <- rep('P',n)
  prob$l  <- rep(0,n)
  
  pred <- prob[,1:numcont]
  trueparm <-matrix(truepar,numcont,1)
  
  meanresp <- trueint+as.matrix(pred)%*%trueparm
  prob$resp<-0
  for (l1 in 1:n)prob$resp[l1] <- rnorm(1,meanresp[l1],respstdev)
  plot(meanresp,prob$resp)
  
  datnow = as.data.frame(prob[,1:numcont])
  fitprob <-lm(prob$resp ~ ., data=datnow)
  
  
  
  ## non prob sample from target
  
  nprob1 <- mvrnorm(n = n1, mu, covm, tol = 1e-6, empirical = FALSE )
  #colnames(nprob1) <- c("X1", "X2")
  nprob1 <- as.data.frame(nprob1)
  nprob1$lab <- rep('NP1',n1)
  nprob1$l  <- rep(1,n1)
  
  pred1 <- nprob1[,1:numcont]
  trueparm <-matrix(truepar,numcont,1)
  
  meanresp1 <- trueint+as.matrix(pred1)%*%trueparm
  nprob1$resp<-0
  for (l1 in 1:n1)nprob1$resp[l1] <- rnorm(1,meanresp1[l1],respstdev)
  
  
  ## non prob sample from poluting
  
  nprob2 <- mvrnorm(n = n2, mu2, covm2, tol = 1e-6, empirical = FALSE )
  #colnames(nprob2) <- c("X1", "X2")
  nprob2 <- as.data.frame(nprob2)
  nprob2$lab <- rep('NP2',n2)
  nprob2$l  <- rep(2,n2)
  
  pred2 <- nprob2[,1:numcont]
  
  meanpoll <- pollint+as.matrix(pred2)%*%pollpar
  nprob2$resp<-0
  for (l1 in 1:n2)nprob2$resp[l1] <- rnorm(1,meanpoll[l1],respstdevpol)
  #plot(meanpoll,nprob2$resp)
  
  
  ### together
  samp <- rbind(prob,nprob1,nprob2)
  
  #scatterplot(samp[,2]~ samp[,1]| lab, data = samp,regLine=FALSE,smooth=FALSE, grid = FALSE, frame = TRUE,
  #            cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="",cex=1.2,pch =c(1,2,19), col = c(1,4,2))
  
  
  new <- rbind(nprob1,nprob2)
  newused <- new
  probused <- prob
  
  
  ### compute  criterion beta change, computes crit as  quantile quantbetachange 
  
  datnoww = as.data.frame(prob[,1:numcont])
  fitprob <-lm(prob$resp ~ ., data=datnoww)
  
  dim <- dim(prob)[1]
  criterion <- matrix(0,dim,1)
  
  for (i in 1:dim) {probnow<- prob[-c(i),]
  datnow = as.data.frame(probnow[,1:numcont])
  fitprobnow <-lm(probnow$resp ~ ., data=datnow)
  #di <- dim(probnow)[1]
  #criterion[i,1] <- abs(fitprobnow$coefficients[2]- fitprob$coefficients[2])/abs(fitprob$coefficients[2]) 
  criterion[i,1] <- Norm(fitprobnow$coefficients- fitprob$coefficients)/Norm(fitprob$coefficients)
  #if ((rs< threshold) & (crit <= critbeta))   {newused$newl[i] <-0}
  } # end i 
  
  #mean(criterion)
  #plot(density(criterion))
  critbeta <- quantile(criterion, quantbetachange)
  
  ####  end compute criterion beta change
  
  
  
    dim <- dim(newused)[1]
    newused$newl <- newused$l
    
    ### loop data +1    
        
    for (i in 1:dim) {columns <- numcont+3
                      probnow<- rbind(prob[,1:columns],newused[i,1:columns])
                      datnow = as.data.frame(probnow[,1:numcont])
                      fitprobnow <-lm(probnow$resp ~ ., data=datnow)
                      #fitprobnow <-lm(X2 ~ X1, data=probnow[,1:2])
                      
                      di <- dim(probnow)[1]
                      rs <- abs(rstudent(fitprobnow)[di])
                      #averagehat <-(sum(hatvalues(fitprobnow))- hatvalues(fitprobnow)[di])/(dim-1)       
                      #maxc <- max(hatvalues(fitprobnow)[1:di])
                      #hatobs <- hatvalues(fitprobnow)[di]
                      #plot(probnow$X1,hatvalues(fitprobnow))
                      crit <- Norm(fitprobnow$coefficients- fitprob$coefficients)/Norm(fitprob$coefficients) 
                      if ((rs< threshold) & (crit <= critbeta)) {newused$newl[i] <-0 }
                      } # end i 
    
    table <- table(newused$l,newused$newl) ### shows how many observations have changed 
    rates[l,1] <-  table[1,1]/sum(table[1,]) ##correctly identified - coming from target
    rates[l,2] <-  table[2,1]/sum(table[2,])
    
    addsample<- newused[ which(newused$newl=='0'),] 
    columns <- numcont+3
    unisamp <- rbind(probused[,1:columns],addsample[,1:columns])
    
  
  ### fit prob
    datnow = as.data.frame(prob[,1:numcont])
    problm <-lm(prob$resp ~ ., data=datnow)
    
    residuals(problm)
  
  ##  fit new
    dats = as.data.frame(unisamp[,1:numcont])
    extlm <-lm(unisamp$resp ~ ., data=dats)
  
  ###
  #results[l,] <- c(problm$coefficients[2],extlm$coefficients[2]) ## estimates slopes prob sample, extended sample
  resultprob[l,] <- problm$coefficients
  resultext[l,] <- extlm$coefficients
  }
####  end loop numsim (data sets)


dev <- matrix(0,numsim,numcont+1)
for (i in 1:numsim)dev[i,]<- resultprob[i,]-truevect
m1<-0
for (i in 1:numsim)m1 <- m1 +Norm(dev[i,])/numsim

devext <- matrix(0,numsim,numcont+1)
for (i in 1:numsim)devext[i,]<- resultext[i,]-truevect
m2<-0
for (i in 1:numsim)m2 <- m2 +Norm(devext[i,])/numsim

resultspss[lp,1]<- m1
resultspss[lp,2]<- m2

} #end probssloop, numsimpss

################################################
#############################################################




############################################################################
# just to look at different sample sizes
plot(numobs,resultspss[,1],ylim=c(0.1,0.6),type='b',lty=3,cex.axis=1.5,cex.lab=1.5,
     main="Norm difference estimated - true for prob sample and extended sample", xlab="Sample size prob sample",
     ylab="")
lines(numobs,resultspss[,2],type='b',lwd=2)
#######################################################



 



#### visualization 

### estimates for all slopes: several boxplots
lim <-numcont+1
for (i in 1:lim){
par <-i
labpar<-c("0", "1","2","3","4")
actlab <- c("Variable",labpar[par])
boxplot(cbind(resultprob[,par],resultext[,par]),names=c("Prob Sample", "Regression Selection"),cex.axis=1.5,
        cex.lab=1.5,cex.main=1.5, main =actlab,cex=1.2)
x<- c(0,4)
y<- c(truevect[par],truevect[par])
lines (x,y)
}

### mean sqrd errors
dev <- matrix(0,numsim,numcont+1)
for (i in 1:numsim)dev[i,]<- resultprob[i,]-truevect
m1<-0
for (i in 1:numsim)m1 <- m1 +Norm(dev[i,])/numsim

devext <- matrix(0,numsim,numcont+1)
for (i in 1:numsim)devext[i,]<- resultext[i,]-truevect
m2<-0
for (i in 1:numsim)m2 <- m2 +Norm(devext[i,])/numsim



cat ("norm err prob sample ", m1)
cat ("norm err extended sample", m2)




#### data last simulation only for one predictor
scatterplot(samp[,4]~ samp[,1]| lab, data = samp,regLine=FALSE,smooth=FALSE,grid = FALSE, frame = TRUE,
            cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Data - only for one predictor meaningful",
            cex=1.2,pch =c(1,2,19), col = c(1,4,2))

trueregr <- trueint + trueslope*samp[,1]
lines(samp[,1],trueregr) 

table(newused$l,newused$newl)


###### new scatter
pred <- samp[,1:numcont]
trueparm <-matrix(truepar,numcont,1)

linpredtrue <- trueint+as.matrix(pred)%*%truepar
scatterplot(samp$resp~ linpredtrue| lab, data = samp,regLine=FALSE,smooth=FALSE,grid = FALSE, frame = TRUE,
            cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Data versus True Predictor",cex=1.2,pch =c(1,2,19), col = c(1,4,2))





##### relative norm
relnorm <- matrix(0,numsim,1)
for (l in 1:numsim) relnorm[l,1]<- Norm(resultprob[l,]-truevect)/Norm(resultext[l,]-truevect)


boxplot(relnorm,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="norm error prob sample / norm error extended sample",
             cex=1.2) #,ylim=c(0, 10))
x<- c(0,4)
y<- c(1,1)
lines (x,y)







### correctly identified
plot(rates[,1],rates[,2],cex.axis=1.5,cex.lab=1.5,cex.main=1.5, 
     main ="Rates of iincluded observations across simulations"
     ,cex=1.2,xlab="Correctly identified as member of  target population",
     ylab="Incorrectly included ")



