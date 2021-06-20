
### regression approach

library(MASS)
library(car)
###propsample

n <- 40 ## probsample
n1 <- 200  ## nonprob from target
n2 <- 400  ## nonprob from poluting

##covariance and mean probsample
var01 <- 1  ## first var
var02 <- 1  ## second var
cov0 <- 0.5
mu01 <- 1
mu02 <- 1
mu <- c(mu01,mu02)
################

##covariance and mean pollution
var21 <- 1  ## first var
var22 <- 1  ## second var
cov2 <-  0.3

mu2 <- c(3,1)
#######################


sigm1<- c(var01,cov0)
sigm2<- c(cov0,var02)
Sigma <- rbind(sigm1,sigm2)



#########################################
#### loop 
numsim <- 30
without <- matrix (0,numsim,1)
with <- matrix (0,numsim,1)
results <- matrix (0,numsim,2)

repet <- 1 #number iterations within
change <- matrix(0,numsim,repet )

threshold <- 1.96 ## for residuals
critbeta <- 0.1   ## for relative change of beta, now computed, replaced by quantbetachange
quantbetachange <- 0.95

for (l in 1: numsim) {
  
  ### prob sample  
  
  prob <- mvrnorm(n , mu, Sigma, tol = 1e-6, empirical = FALSE )
  colnames(prob) <- c("X1", "X2")
  prob <- as.data.frame(prob)
  prob$lab <- rep('P',n)
  prob$l  <- rep(0,n)
  
  
  ## non prob sample from target
  
  nprob1 <- mvrnorm(n = n1, mu, Sigma, tol = 1e-6, empirical = FALSE )
  colnames(nprob1) <- c("X1", "X2")
  nprob1 <- as.data.frame(nprob1)
  nprob1$lab <- rep('NP1',n1)
  nprob1$l  <- rep(1,n1)
  
  
  ## non prob sample from poluting
  sigm1<- c(var21,cov2)
  sigm2<- c(cov2,var22)
  Sigma2 <- rbind(sigm1,sigm2)
  nprob2 <- mvrnorm(n = n2, mu2, Sigma2, tol = 1e-6, empirical = FALSE )
  colnames(nprob2) <- c("X1", "X2")
  nprob2 <- as.data.frame(nprob2)
  nprob2$lab <- rep('NP2',n2)
  nprob2$l  <- rep(2,n2)
  
  
  
  ### together
  samp <- rbind(prob,nprob1,nprob2)
  
  #scatterplot(samp[,2]~ samp[,1]| lab, data = samp,regLine=FALSE,smooth=FALSE, grid = FALSE, frame = TRUE,
  #            cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="",cex=1.2,pch =c(1,2,19), col = c(1,4,2))
  
  
  new <- rbind(nprob1,nprob2)
  newused <- new
  probused <- prob
  
  
  ### compute  criterion beta change, computes crit as .95 quantile
  fitprob <-lm(X2 ~ X1, data=prob[,1:2])
  dim <- dim(prob)[1]
  
  criterion <- matrix(0,dim,1)
  for (i in 1:dim) {probnow<- prob[-c(i),]
  fitprobnow <-lm(X2 ~ X1, data=probnow[,1:2])
  di <- dim(probnow)[1]
  criterion[i,1] <- abs(fitprobnow$coefficients[2]- fitprob$coefficients[2])/abs(fitprob$coefficients[2]) 
  if ((rs< threshold) & (crit <= critbeta)) 
    #if (rs< threshold)   
  {newused$newl[i] <-0 
  }} # end i 
  #mean(criterion)
  #plot(density(criterion))
  critbeta <- quantile(criterion, quantbetachange)
  ####  end compute criterion beta change
  
  
  #### additional loop prob, now omitted
  
  #for( add in 1:repet) { 
    
    #m0 <- mahalanobis(probused[,1:2],colMeans(probused[,1:2]),cov(probused[,1:2])) 
    #plot (density(m0))  
    #q <- quantile(m0, 0.7)
    
    #mnew <- mahalanobis(newused[,1:2],colMeans(probused[,1:2]),cov(probused[,1:2])) 
    #plot (density(mnew))
    fitprob <-lm(X2 ~ X1, data=prob[,1:2])
    dim <- dim(newused)[1]
    newused$newl <- newused$l
    
    ### loop data +1    
        
    for (i in 1:dim) {probnow<- rbind(prob[,1:2],newused[i,1:2])
                      fitprobnow <-lm(X2 ~ X1, data=probnow[,1:2])
                      di <- dim(probnow)[1]
                      rs <- abs(rstudent(fitprobnow)[di])
                      #averagehat <-(sum(hatvalues(fitprobnow))- hatvalues(fitprobnow)[di])/(dim-1)       
                      #maxc <- max(hatvalues(fitprobnow)[1:di])
                      #hatobs <- hatvalues(fitprobnow)[di]
                      #plot(probnow$X1,hatvalues(fitprobnow))
                      crit <- abs(fitprobnow$coefficients[2]- fitprob$coefficients[2])/abs(fitprob$coefficients[2]) 
                      if ((rs< threshold) & (crit <= critbeta)) 
                        #if (rs< threshold)   
                        {newused$newl[i] <-0 
            }} # end i 
    
    table(newused$l,newused$newl) ### shows how many observations have changed 
    
    addsample<- newused[ which(newused$newl=='0'),] 
    unisamp <- rbind(probused[,1:2],addsample[,1:2])
    
    #probused <- unisamp
    #newused <- newused[ which(newused$newl!='0'),]
    
 # } #add end
  
  
  
  ####
  
  
  
  ### fit prob
  problm <-lm(X2 ~ X1, data=prob[,1:2])
  residuals(problm)
  
  ##  fit new
  
  extlm <-lm(X2 ~ X1, data=unisamp)
  ###
  results[l,] <- c(problm$coefficients[2],extlm$coefficients[2])
}
####  end loop


#### visualization 

boxplot(results,names=c("Prob Sample", "Regression Selection"),cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="",cex=1.2)
x<- c(0,4)
y<- c(trueslope,trueslope)
lines (x,y)

mean((results[,1]-trueslope)^2)
mean((results[,2]-trueslope)^2)

mean(abs(results[,1]-trueslope))
mean(abs(results[,2]-trueslope))

scatterplot(samp[,2]~ samp[,1]| lab, data = samp,regLine=FALSE,smooth=FALSE,
            grid = FALSE, frame = TRUE,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="",cex=1.2,pch =c(1,2,19),
            col = c(1,4,2))


table(newused$l,newused$newl)

rel <- (results[,1]-trueslope)^2/((results[,2]-trueslope)^2)
boxplot(rel,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="sqrt error prob sample / error extended sample",cex=1.2
        ,ylim=c(0, 10))
x<- c(0,4)
y<- c(1,1)
lines (x,y)



rel2 <- abs(results[,1]-trueslope)/(abs(results[,2]-trueslope))
boxplot(rel2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="abs error prob sample / error extended sample",cex=1.2,
        ylim=c(0, 5))
x<- c(0,4)
y<- c(1,1)
lines (x,y)



#### trial: select threeshold for beta change now included in  program above

fitprob <-lm(X2 ~ X1, data=prob[,1:2])
dim <- dim(prob)[1]

criterion <- matrix(0,dim,1)
for (i in 1:dim) {probnow<- prob[-c(i),]
fitprobnow <-lm(X2 ~ X1, data=probnow[,1:2])
di <- dim(probnow)[1]
criterion[i,1] <- abs(fitprobnow$coefficients[2]- fitprob$coefficients[2])/abs(fitprob$coefficients[2]) 
if ((rs< threshold) & (crit <= critbeta)) 
  #if (rs< threshold)   
{newused$newl[i] <-0 
}} # end i 
mean(criterion)
plot(density(criterion))
quantile(criterion, 0.95)

