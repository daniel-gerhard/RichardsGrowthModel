#############################
library(truncnorm)

richards <- function(x, beta){
  mu <- beta[1] * (1 + (beta[2]-1) * exp(-beta[3]*(x-beta[4])))^(1/(1-beta[2]))
  mu[is.na(mu)] <- 0
  return(mu)
}

############################################################
## MH
richardsMH1 <- function(x, y, bmu, bsd, atau, btau, betastart, iter){
  loglik <- function(x, y, beta){
    mu <- richards(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(beta[5]), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- richards(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  likp <- function(x, y, beta){
    mu <- richards(x, beta)
    dnorm(y, mu, 1/sqrt(beta[5]))
  }
  
  pm <- matrix(NA, ncol=length(betastart), nrow=iter)
  pm[1,] <- betastart
  accept <- numeric(length=4)
  
  llk <- numeric(length=iter)
  llk[1] <- loglik(x, y, betastart)
  
  lmat <- matrix(NA, ncol=length(y), nrow=iter)
  lmat[1,] <- likp(x, y, betastart)
  
  for (i in 2:iter){
    bn <- bo <- pm[i-1, ]
    
    ### sampling a
    bn[1] <- rtruncnorm(1, a=0, mean=bo[1], sd=1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[1], a=0, mean=bmu[1], sd=bsd[1])) - 
      log(dtruncnorm(bo[1], a=0, mean=bmu[1], sd=bsd[1])) -
      log(dtruncnorm(bn[1], a=0, mean=bo[1], sd=1)) +
      log(dtruncnorm(bo[1], a=0, mean=bn[1], sd=1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[1] <- accept[1] + 1
    } else {
      bn <- bo
    }
    
    ### sampling d1
    bn[2] <- rtruncnorm(1, a=1, mean=bo[2], sd=0.5)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[2], a=1, mean=bmu[2], sd=bsd[2])) - 
      log(dtruncnorm(bo[2], a=1, mean=bmu[2], sd=bsd[2])) -
      log(dtruncnorm(bn[2], a=1, mean=bo[2], sd=0.5)) +
      log(dtruncnorm(bo[2], a=1, mean=bn[2], sd=0.5))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[2] <- accept[2] + 1
    } else {
      bn <- bo
    }
    
    ### sampling k1
    bn[3] <- rtruncnorm(1, a=0, mean=bo[3], sd=0.5)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[3], a=0, mean=bmu[3], sd=bsd[3])) - 
      log(dtruncnorm(bo[3], a=0, mean=bmu[3], sd=bsd[3])) -
      log(dtruncnorm(bn[3], a=0, mean=bo[3], sd=0.5)) +
      log(dtruncnorm(bo[3], a=0, mean=bn[3], sd=0.5))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[3] <- accept[3] + 1
    } else {
      bn <- bo
    }
    
    ### sampling g1
    bn[4] <- rtruncnorm(1, a=0, mean=bo[4], sd=1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[4], a=0, mean=bmu[4], sd=bsd[4])) - 
      log(dtruncnorm(bo[4], a=0, mean=bmu[4], sd=bsd[4])) -
      log(dtruncnorm(bn[4], a=0, mean=bo[4], sd=1)) +
      log(dtruncnorm(bo[4], a=0, mean=bn[4], sd=1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[4] <- accept[4] + 1
    } else {
      bn <- bo
    } 
    
    ### sampling tau (Gibbs)
    bn[5] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
    pm[i,] <- bn
    
    ### likelihood
    lmat[i,] <- likp(x, y, bn)
  }
  paccept <- accept/iter
  
  res <- list(pm=pm, lmat=lmat, paccept=paccept)
  return(res)
}

##################

richardsMH2 <- function(x, y, bmu, bsd, atau, btau, betastart, iter){
  loglik <- function(x, y, beta){
    mu <- richards(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(beta[5]), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- richards(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  likp <- function(x, y, beta){
    mu <- richards(x, beta)
    dnorm(y, mu, 1/sqrt(beta[5]))
  }
  
  pm <- matrix(NA, ncol=length(betastart), nrow=iter)
  pm[1,] <- betastart
  accept <- numeric(length=4)
  
  llk <- numeric(length=iter)
  llk[1] <- loglik(x, y, betastart)
  
  lmat <- matrix(NA, ncol=length(y), nrow=iter)
  lmat[1,] <- likp(x, y, betastart)
  
  for (i in 2:iter){
    bn <- bo <- pm[i-1, ]
    
    ### sampling a
    bn[1] <- rtruncnorm(1, a=0, mean=bo[1], sd=1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[1], a=0, mean=bmu[1], sd=bsd[1])) - 
      log(dtruncnorm(bo[1], a=0, mean=bmu[1], sd=bsd[1])) -
      log(dtruncnorm(bn[1], a=0, mean=bo[1], sd=1)) +
      log(dtruncnorm(bo[1], a=0, mean=bn[1], sd=1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[1] <- accept[1] + 1
    } else {
      bn <- bo
    }
    
    ### sampling d2
    bn[2] <- rtruncnorm(1, a=-1, b=1, mean=bo[2], sd=0.1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[2], a=-1, b=1, mean=bmu[2], sd=bsd[2])) - 
      log(dtruncnorm(bo[2], a=-1, b=1, mean=bmu[2], sd=bsd[2])) -
      log(dtruncnorm(bn[2], a=-1, b=1, mean=bo[2], sd=0.1)) +
      log(dtruncnorm(bo[2], a=-1, b=1, mean=bn[2], sd=0.1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[2] <- accept[2] + 1
    } else {
      bn <- bo
    }
    
    ### sampling k2
    bn[3] <- rtruncnorm(1, a=0, mean=bo[3], sd=0.5)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[3], a=0, mean=bmu[3], sd=bsd[3])) - 
      log(dtruncnorm(bo[3], a=0, mean=bmu[3], sd=bsd[3])) -
      log(dtruncnorm(bn[3], a=0, mean=bo[3], sd=0.5)) +
      log(dtruncnorm(bo[3], a=0, mean=bn[3], sd=0.5))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[3] <- accept[3] + 1
    } else {
      bn <- bo
    }  
    
    ### sampling g2
    bn[4] <- rtruncnorm(1, a=0, mean=bo[4], sd=1)
    logR <- loglik(x, y, bn) -
      loglik(x, y, bo) +
      log(dtruncnorm(bn[4], a=0, mean=bmu[4], sd=bsd[4])) -
      log(dtruncnorm(bo[4], a=0, mean=bmu[4], sd=bsd[4])) -
      log(dtruncnorm(bn[4], a=0, mean=bo[4], sd=1)) +
      log(dtruncnorm(bo[4], a=0, mean=bn[4], sd=1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[4] <- accept[4] + 1
    } else {
      bn <- bo
    }
    
    ### sampling tau (Gibbs)
    bn[5] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
    pm[i,] <- bn
    
    ### likelihood
    lmat[i,] <- likp(x, y, bn)
  }
  paccept <- accept/iter
  
  res <- list(pm=pm, lmat=lmat, paccept=paccept)
  return(res)
}

###############3

richardsMH3 <- function(x, y, bmu, bsd, atau, btau, betastart, iter){
  loglik <- function(x, y, beta){
    mu <- richards(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(beta[5]), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- richards(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  likp <- function(x, y, beta){
    mu <- richards(x, beta)
    dnorm(y, mu, 1/sqrt(beta[5]))
  }
  
  pm <- matrix(NA, ncol=length(betastart), nrow=iter)
  pm[1,] <- betastart
  accept <- numeric(length=4)
  
  llk <- numeric(length=iter)
  llk[1] <- loglik(x, y, betastart)
  
  lmat <- matrix(NA, ncol=length(y), nrow=iter)
  lmat[1,] <- likp(x, y, betastart)
  
  for (i in 2:iter){
    bn <- bo <- pm[i-1, ]
    
    ### sampling a
    bn[1] <- rtruncnorm(1, a=0, mean=bo[1], sd=1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[1], a=0, mean=bmu[1], sd=bsd[1])) - 
      log(dtruncnorm(bo[1], a=0, mean=bmu[1], sd=bsd[1])) -
      log(dtruncnorm(bn[1], a=0, mean=bo[1], sd=0.5)) +
      log(dtruncnorm(bo[1], a=0, mean=bn[1], sd=0.5))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[1] <- accept[1] + 1
    } else {
      bn <- bo
    }
    
    
    ### sampling k1
    bn[3] <- rtruncnorm(1, a=0, mean=bo[3], sd=0.5)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[3], a=0, mean=bmu[3], sd=bsd[3])) - 
      log(dtruncnorm(bo[3], a=0, mean=bmu[3], sd=bsd[3])) -
      log(dtruncnorm(bn[3], a=0, mean=bo[3], sd=0.5)) +
      log(dtruncnorm(bo[3], a=0, mean=bn[3], sd=0.5))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[3] <- accept[3] + 1
    } else {
      bn <- bo
    }
    
    ### sampling g1
    bn[4] <- rtruncnorm(1, a=0, mean=bo[4], sd=1)
    logR <- loglik(x, y, bn) - 
      loglik(x, y, bo) +
      log(dtruncnorm(bn[4], a=0, mean=bmu[4], sd=bsd[4])) - 
      log(dtruncnorm(bo[4], a=0, mean=bmu[4], sd=bsd[4])) -
      log(dtruncnorm(bn[4], a=0, mean=bo[4], sd=1)) +
      log(dtruncnorm(bo[4], a=0, mean=bn[4], sd=1))
    logU <- log(runif(1, 0, 1))
    if (logR > logU){
      bo <- bn
      accept[4] <- accept[4] + 1
    } else {
      bn <- bo
    } 
    
    ### sampling tau (Gibbs)
    bn[5] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
    pm[i,] <- bn
    
    ### likelihood
    lmat[i,] <- likp(x, y, bn)
  }
  paccept <- accept/iter
  
  res <- list(pm=pm, lmat=lmat, paccept=paccept)
  return(res)
}