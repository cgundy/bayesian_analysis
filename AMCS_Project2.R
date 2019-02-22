
# using Bayesian inference to predict likelihood of contraceptive use amond Indonesian women

rm(list=ls())
require(mvtnorm)
library(MASS)
cmc <- read.table("ACMS_Project2.txt", sep = ",", header = FALSE)
names(cmc) <- c("wife_age","wife_ed","husband_ed","children","wife_rel",
                "wife_work","husband_oc","SOV","media","contraceptive")
cmc$contraceptive <- as.factor(cmc$contraceptive)
levels(cmc$contraceptive) <- c(0,1,1)
cmc$contraceptive <- as.numeric(cmc$contraceptive)
cmc$contraceptive <- (cmc$contraceptive-1)
set.seed(3)
cmc<- cmc[sample(nrow(cmc),900),]

#obtain MLEs for grade data
attach(cmc)
mod.mle <- glm(contraceptive~wife_age+wife_ed+children+wife_rel+wife_work,family=binomial)
set.seed(3)


#Write function to compute likelihoods. Input betas and data.
likelihood <- function(betas,y,x1,x2,x3,x4,x5){
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  beta4 <- betas[5]
  beta5 <- betas[6]
  #compute fitted values
  p.hats <- plogis(beta0+beta1*x1+beta2*x2+beta3*x3+beta4*x4+beta5*x5)
  #compute log-likelihood 
  #log.like <- sum(y*log(p.hats)+(1-y)*log(1-p.hats))
  #compute likelihood 
  like <- prod(p.hats^y*(1-p.hats)^(1-y))
}


#Set parameters of prior distribution
prior.mean.b0 <- 0
prior.mean.b1 <- 0
prior.mean.b2 <- 0
prior.mean.b3 <- 0
prior.mean.b4 <- 0
prior.mean.b5 <- 0
prior.var.b0 <- 10000
prior.var.b1 <- 10000
prior.var.b2 <- 10000
prior.var.b3 <- 10000
prior.var.b4 <- 10000
prior.var.b5 <- 10000
prior.mu <- c(prior.mean.b0,prior.mean.b1,prior.mean.b2,prior.mean.b3,prior.mean.b4,prior.mean.b5)
prior.sigma <- diag(c(prior.var.b0,prior.var.b1,prior.var.b2,prior.var.b3,prior.var.b4,prior.var.b5))


#STEP 1: pick starting values, i.e., use MLEs from standard model as starting values
betas.current <- coef(mod.mle)

#Count number of times proposed point is accepted
accept <- 0
#Set up NULL vectors in which to store posterior samples
posterior.b0 <- NULL
posterior.b1 <- NULL
posterior.b2 <- NULL
posterior.b3 <- NULL
posterior.b4 <- NULL
posterior.b5 <- NULL
#MH parameter to adjust acceptance rate
sigma.MH <- 0.75
for (j in 1:5000){
  #STEP 2: Use current betas to compute posterior (with flat prior)
  post.current <- likelihood(betas.current,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*
    dmvnorm(betas.current,prior.mu,prior.sigma)
  #STEP 3: Generate new propsed betas from multivariate normal
  betas.proposed <- mvrnorm(n = 1, mu=c(betas.current[1],betas.current[2],betas.current[3],betas.current[4],betas.current[5],betas.current[6]), Sigma=sigma.MH*vcov(mod.mle))
  #STEP 4: Compute posterior with proposed betas
  post.proposed <- likelihood(betas.proposed,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*
    dmvnorm(betas.proposed,prior.mu,prior.sigma)
  #STEP 5: Compute probability of moving
  pmove <- min(post.proposed/post.current,1)
  #STEP 6: Generate random uniform(0,1) 	    
  u <- runif(1,0,1)
  #STEP 7: Move to proposed parameter values if u<pmove
  if (u < pmove){
    betas.current <- betas.proposed
    accept <- accept + 1
  }
  posterior.b0 <- c(posterior.b0,betas.current[1])
  posterior.b1 <- c(posterior.b1,betas.current[2])
  posterior.b2 <- c(posterior.b2,betas.current[3])
  posterior.b3 <- c(posterior.b3,betas.current[4])
  posterior.b4 <- c(posterior.b4,betas.current[5])
  posterior.b5 <- c(posterior.b5,betas.current[6])
  
}

accept/5000

summary(posterior.b0[1001:5000])
summary(posterior.b1[1001:5000])
summary(posterior.b2[1001:5000])
summary(posterior.b3[1001:5000])
summary(posterior.b4[1001:5000])
summary(posterior.b5[1001:5000])

##Do early values represent the true posterior? 
#Constuct trace plot as in Figure 10.1 of text
plot(1:5000,posterior.b1,type="l",xlab="Iteration",ylab="Beta1 Values")
par(mfrow=c(2,2))
plot(1:5000,posterior.b2,type="l",xlab="Iteration",ylab="Beta2 Values")
plot(1:5000,posterior.b3,type="l",xlab="Iteration",ylab="Beta3 Values")
plot(1:5000,posterior.b4,type="l",xlab="Iteration",ylab="Beta4 Values")
plot(1:5000,posterior.b5,type="l",xlab="Iteration",ylab="Beta5 Values")

par(mfrow=c(1,1))

prior.mean.b0 <- 0
prior.mean.b1 <- -1
prior.mean.b2 <- 1
prior.mean.b3 <- 0
prior.mean.b4 <- -1
prior.mean.b5 <- 1
prior.var.b0 <- 100
prior.var.b1 <- 100
prior.var.b2 <- 10
prior.var.b3 <- 100
prior.var.b4 <- 10
prior.var.b5 <- 100
prior.mu <- c(prior.mean.b0,prior.mean.b1,prior.mean.b2,prior.mean.b3,prior.mean.b4,prior.mean.b5)
prior.sigma <- diag(c(prior.var.b0,prior.var.b1,prior.var.b2,prior.var.b3,prior.var.b4,prior.var.b5))



#STEP 1: pick starting values for EACH chain; vary the slopes a bit
betas.current1 <- coef(mod.mle)
betas.current2 <- c(coef(mod.mle)[1],coef(mod.mle)[2]+0.03,coef(mod.mle)[3]+0.03,coef(mod.mle)[4]+0.03,coef(mod.mle)[5]+0.03,coef(mod.mle)[6]+0.03)
betas.current3 <- c(coef(mod.mle)[1],coef(mod.mle)[2]-0.03,coef(mod.mle)[3]-0.03,coef(mod.mle)[4]-0.03,coef(mod.mle)[5]-0.03,coef(mod.mle)[6]-0.03)

#Count number of times proposed point is accepted
accept1 <- 0
accept2 <- 0
accept3 <- 0
#Set up NULL vectors in which to store posterior samples
posterior1.b0 <- NULL
posterior1.b1 <- NULL
posterior1.b2 <- NULL
posterior1.b3 <- NULL
posterior1.b4 <- NULL
posterior1.b5 <- NULL
posterior2.b0 <- NULL
posterior2.b1 <- NULL
posterior2.b2 <- NULL
posterior2.b3 <- NULL
posterior2.b4 <- NULL
posterior2.b5 <- NULL
posterior3.b0 <- NULL
posterior3.b1 <- NULL
posterior3.b2 <- NULL
posterior3.b3 <- NULL
posterior3.b4 <- NULL
posterior3.b5 <- NULL
#MH parameter to adjust acceptance rate
sigma.MH1 <- 0.75
sigma.MH2 <- 0.75
sigma.MH3 <- 0.75

for (j in 1:5000){
  post.current1 <- likelihood(betas.current1,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.current1,prior.mu,prior.sigma)
  post.current2 <- likelihood(betas.current2,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.current2,prior.mu,prior.sigma)
  post.current3 <- likelihood(betas.current3,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.current3,prior.mu,prior.sigma)
  #STEP 3: Generate new propsed betas from multivariate normal
  betas.proposed1 <- mvrnorm(n = 1, mu=c(betas.current1[1],betas.current1[2],betas.current1[3],betas.current1[4],betas.current1[5],betas.current1[6]), Sigma=sigma.MH*vcov(mod.mle))
  betas.proposed2 <- mvrnorm(n = 1, mu=c(betas.current2[1],betas.current2[2],betas.current2[3],betas.current2[4],betas.current2[5],betas.current2[6]), Sigma=sigma.MH*vcov(mod.mle))
  betas.proposed3 <- mvrnorm(n = 1, mu=c(betas.current3[1],betas.current3[2],betas.current3[3],betas.current3[4],betas.current3[5],betas.current3[6]), Sigma=sigma.MH*vcov(mod.mle))
  #STEP 4: Compute posterior with proposed betas
  post.proposed1 <- likelihood(betas.proposed1,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.proposed1,prior.mu,prior.sigma)
  post.proposed2 <- likelihood(betas.proposed2,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.proposed2,prior.mu,prior.sigma)
  post.proposed3 <- likelihood(betas.proposed3,contraceptive,wife_age,wife_ed,children,wife_rel,wife_work)*dmvnorm(betas.proposed3,prior.mu,prior.sigma)
  #STEP 5: Compute probability of moving
  pmove1 <- min(post.proposed1/post.current1,1)
  pmove2 <- min(post.proposed2/post.current2,1)
  pmove3 <- min(post.proposed3/post.current3,1)
  #STEP 6: Generate random uniform(0,1) 	    
  u1 <- runif(1,0,1)
  u2 <- runif(1,0,1)
  u3 <- runif(1,0,1)
  #STEP 7: Move to proposed parameter values if u<pmove
  if (u1 < pmove1){
    betas.current1 <- betas.proposed1
    accept1 <- accept1 + 1
  }
  if (u2 < pmove2){
    betas.current2 <- betas.proposed2
    accept2 <- accept2 + 1
  }
  if (u3 < pmove3){
    betas.current3 <- betas.proposed3
    accept3 <- accept3 + 1
  }
  #Chain 1
  posterior1.b0 <- c(posterior1.b0,betas.current1[1])
  posterior1.b1 <- c(posterior1.b1,betas.current1[2])
  posterior1.b2 <- c(posterior1.b2,betas.current1[3])
  posterior1.b3 <- c(posterior1.b3,betas.current1[4])
  posterior1.b4 <- c(posterior1.b4,betas.current1[5])
  posterior1.b5 <- c(posterior1.b5,betas.current1[6])
  #Chain 2
  posterior2.b0 <- c(posterior2.b0,betas.current2[1])
  posterior2.b1 <- c(posterior2.b1,betas.current2[2])
  posterior2.b2 <- c(posterior2.b2,betas.current2[3])
  posterior2.b3 <- c(posterior2.b3,betas.current2[4])
  posterior2.b4 <- c(posterior2.b4,betas.current2[5])
  posterior2.b5 <- c(posterior2.b5,betas.current2[6])
  #Chain 3
  posterior3.b0 <- c(posterior3.b0,betas.current3[1])
  posterior3.b1 <- c(posterior3.b1,betas.current3[2])
  posterior3.b2 <- c(posterior3.b2,betas.current3[3])
  posterior3.b3 <- c(posterior3.b3,betas.current3[4])
  posterior3.b4 <- c(posterior3.b4,betas.current3[5])
  posterior3.b5 <- c(posterior3.b5,betas.current3[6])
}

accept1/5000
accept2/5000
accept3/5000

#Constuct trace plot for all three chains as in Figure 10.3 of text
plot(1:5000,posterior1.b1,type="l",xlab="Iteration",ylab="Simulated Beta1 Values",col="red")
lines(posterior2.b1,type="l",col="blue")
lines(posterior3.b1,type="l",col="green")

plot(1:5000,posterior1.b2,type="l",xlab="Iteration",ylab="Simulated Beta2 Values",col="red")
lines(posterior2.b2,type="l",col="blue")
lines(posterior3.b2,type="l",col="green")

plot(1:5000,posterior1.b3,type="l",xlab="Iteration",ylab="Simulated Beta3 Values",col="red")
lines(posterior2.b3,type="l",col="blue")
lines(posterior3.b3,type="l",col="green")

plot(1:5000,posterior1.b4,type="l",xlab="Iteration",ylab="Simulated Beta4 Values",col="red")
lines(posterior2.b4,type="l",col="blue")
lines(posterior3.b4,type="l",col="green")

plot(1:5000,posterior1.b5,type="l",xlab="Iteration",ylab="Simulated Beta5 Values",col="red")
lines(posterior2.b5,type="l",col="blue")
lines(posterior3.b5,type="l",col="green")


###Compute BGR statistic for assessing convergence (pages 417-418)
#compute variance of each chain
var1.b1 <- var(posterior1.b1)
var2.b1 <- var(posterior2.b1)
var3.b1 <- var(posterior3.b1)
var1.b2 <- var(posterior1.b2)
var2.b2 <- var(posterior2.b2)
var3.b2 <- var(posterior3.b2)
var1.b3 <- var(posterior1.b3)
var2.b3 <- var(posterior2.b3)
var3.b3 <- var(posterior3.b3)
var1.b4 <- var(posterior1.b4)
var2.b4 <- var(posterior2.b4)
var3.b4 <- var(posterior3.b4)
var1.b5 <- var(posterior1.b5)
var2.b5 <- var(posterior2.b5)
var3.b5 <- var(posterior3.b5)
#compute W, the average of the variances
W.b1 <- mean(c(var1.b1,var2.b1,var3.b1))
W.b2 <- mean(c(var1.b2,var2.b2,var3.b2))
W.b3 <- mean(c(var1.b3,var2.b3,var3.b3))
W.b4 <- mean(c(var1.b4,var2.b4,var3.b4))
W.b5 <- mean(c(var1.b5,var2.b5,var3.b5))
##compute B, the between chain variability. 
#Compute grand mean, theta-bar-bar and mean of each chain
theta.bar.bar.b1 <- mean(c(posterior1.b1,posterior2.b1,posterior3.b1))
theta.bar.bar.b2 <- mean(c(posterior1.b2,posterior2.b2,posterior3.b2))
theta.bar.bar.b3 <- mean(c(posterior1.b3,posterior2.b3,posterior3.b3))
theta.bar.bar.b4 <- mean(c(posterior1.b4,posterior2.b4,posterior3.b4))
theta.bar.bar.b5 <- mean(c(posterior1.b5,posterior2.b5,posterior3.b5))

mean1.b1 <- mean(posterior1.b1)
mean2.b1 <- mean(posterior2.b1)
mean3.b1 <- mean(posterior3.b1)
mean1.b2 <- mean(posterior1.b2)
mean2.b2 <- mean(posterior2.b2)
mean3.b2 <- mean(posterior3.b2)
mean1.b3 <- mean(posterior1.b3)
mean2.b3 <- mean(posterior2.b3)
mean3.b3 <- mean(posterior3.b3)
mean1.b4 <- mean(posterior1.b4)
mean2.b4 <- mean(posterior2.b4)
mean3.b4 <- mean(posterior3.b4)
mean1.b5 <- mean(posterior1.b5)
mean2.b5 <- mean(posterior2.b5)
mean3.b5 <- mean(posterior3.b5)

#Compute B 
B.b1 <- 0.5*(5000*(mean1.b1-theta.bar.bar.b1)^2+5000*(mean2.b1-theta.bar.bar.b1)^2+5000*(mean3.b1-theta.bar.bar.b1)^2)
B.b2 <- 0.5*(5000*(mean1.b2-theta.bar.bar.b2)^2+5000*(mean2.b2-theta.bar.bar.b2)^2+5000*(mean3.b2-theta.bar.bar.b2)^2)
B.b3 <- 0.5*(5000*(mean1.b3-theta.bar.bar.b3)^2+5000*(mean2.b3-theta.bar.bar.b3)^2+5000*(mean3.b3-theta.bar.bar.b3)^2)
B.b4 <- 0.5*(5000*(mean1.b4-theta.bar.bar.b4)^2+5000*(mean2.b4-theta.bar.bar.b4)^2+5000*(mean3.b4-theta.bar.bar.b4)^2)
B.b5 <- 0.5*(5000*(mean1.b5-theta.bar.bar.b5)^2+5000*(mean2.b5-theta.bar.bar.b5)^2+5000*(mean3.b5-theta.bar.bar.b5)^2)
#Compute V.hat
V.hat.b1 <- (4999/5000)*W.b1+(1/5000)*B.b1
V.hat.b2 <- (4999/5000)*W.b2+(1/5000)*B.b2
V.hat.b3 <- (4999/5000)*W.b3+(1/5000)*B.b3
V.hat.b4 <- (4999/5000)*W.b4+(1/5000)*B.b4
V.hat.b5 <- (4999/5000)*W.b5+(1/5000)*B.b5
#Compute R.hat
R.hat.b1 <- sqrt(V.hat.b1/W.b1)
R.hat.b2 <- sqrt(V.hat.b2/W.b2)
R.hat.b3 <- sqrt(V.hat.b3/W.b3)
R.hat.b4 <- sqrt(V.hat.b4/W.b4)
R.hat.b5 <- sqrt(V.hat.b5/W.b5)

