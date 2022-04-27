library(truncnorm)

# Double Richards function
drichards <- function(x, beta){
  mu1 <- beta[1] * (1 + (beta[2]-1) * exp(-beta[3]*(x-beta[4])))^(1/(1-beta[2]))
  mu1[is.na(mu1)] <- 0
  mu2 <- beta[5] * (1 + (beta[6]-1) * exp(-beta[7]*(x-(beta[4]+beta[8]))))^(1/(1-beta[6]))
  mu2[is.na(mu2)] <- 0
  fm <- mu1 + mu2
  return(fm)
}

# Double Richards Sampler
richardsMH <- function(x, y, bmu, bsd, atau, btau, betastart, iter, ltd1, utd1, ltd2, utd2, fixed){
  loglik <- function(x, y, beta){
    mu <- drichards(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(beta[9]), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- drichards(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  likp <- function(x, y, beta){
    mu <- drichards(x, beta)
    dnorm(y, mu, 1/sqrt(beta[9]))
  }
  pm <- matrix(NA, ncol=length(betastart), nrow=iter)
  pm[1,] <- betastart
  accept <- numeric(length=8)
  
  llk <- numeric(length=iter)
  llk[1] <- loglik(x, y, betastart)
  
  lmat <- matrix(NA, ncol=length(y), nrow=iter)
  lmat[1,] <- likp(x, y, betastart)
  
  for (i in 2:iter){
    bn <- bo <- pm[i-1, ]
    ### sampling a1
    if (is.na(fixed[1])){
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
    }
    ### sampling d1
    if (is.na(fixed[2])){
      bn[2] <- rtruncnorm(1, a=ltd1, b=utd1, mean=bo[2], sd=0.5)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[2], a=ltd1, b=utd1, mean=bmu[2], sd=bsd[2])) - 
        log(dtruncnorm(bo[2], a=ltd1, b=utd1, mean=bmu[2], sd=bsd[2])) -
        log(dtruncnorm(bn[2], a=ltd1, b=utd1, mean=bo[2], sd=0.5)) +
        log(dtruncnorm(bo[2], a=ltd1, b=utd1, mean=bn[2], sd=0.5))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        accept[2] <- accept[2] + 1
      } else {
        bn <- bo
      }
    }
    ### sampling k1
    if (is.na(fixed[3])){
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
    }
    ### sampling g1
    if (is.na(fixed[4])){
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
    }
    ### sampling a2
    if (is.na(fixed[5])){
      bn[5] <- rtruncnorm(1, a=0, mean=bo[5], sd=1)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[5], a=0, mean=bmu[5], sd=bsd[5])) - 
        log(dtruncnorm(bo[5], a=0, mean=bmu[5], sd=bsd[5])) -
        log(dtruncnorm(bn[5], a=0, mean=bo[5], sd=1)) +
        log(dtruncnorm(bo[5], a=0, mean=bn[5], sd=1))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        accept[5] <- accept[5] + 1
      } else {
        bn <- bo
      }
    }
    ### sampling d2
    if (is.na(fixed[6])){
      bn[6] <- rtruncnorm(1, a=ltd2, b=utd2, mean=bo[6], sd=0.5)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[6], a=ltd2, b=utd2, mean=bmu[6], sd=bsd[6])) - 
        log(dtruncnorm(bo[6], a=ltd2, b=utd2, mean=bmu[6], sd=bsd[6])) -
        log(dtruncnorm(bn[6], a=ltd2, b=utd2, mean=bo[6], sd=0.5)) +
        log(dtruncnorm(bo[6], a=ltd2, b=utd2, mean=bn[6], sd=0.5))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        accept[2] <- accept[2] + 1
      } else {
        bn <- bo
      }
    }
    ### sampling k2
    if (is.na(fixed[7])){
      bn[7] <- rtruncnorm(1, a=0, mean=bo[7], sd=0.5)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[7], a=0, mean=bmu[7], sd=bsd[7])) - 
        log(dtruncnorm(bo[7], a=0, mean=bmu[7], sd=bsd[7])) -
        log(dtruncnorm(bn[7], a=0, mean=bo[7], sd=0.5)) +
        log(dtruncnorm(bo[7], a=0, mean=bn[7], sd=0.5))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        accept[7] <- accept[7] + 1
      } else {
        bn <- bo
      }
    }
    ### sampling g2
    if (is.na(fixed[8])){
      bn[8] <- rtruncnorm(1, a=0, mean=bo[8], sd=1)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[8], a=0, mean=bmu[8], sd=bsd[8])) - 
        log(dtruncnorm(bo[8], a=0, mean=bmu[8], sd=bsd[8])) -
        log(dtruncnorm(bn[8], a=0, mean=bo[8], sd=1)) +
        log(dtruncnorm(bo[8], a=0, mean=bn[8], sd=1))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        accept[8] <- accept[8] + 1
      } else {
        bn <- bo
      } 
    }
    ### sampling tau (Gibbs)
    bn[9] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
    pm[i,] <- bn
    
    ### likelihood
    lmat[i,] <- likp(x, y, bn)
  }
  paccept <- accept/iter
  
  res <- list(pm=pm, lmat=lmat, paccept=paccept)
  return(res)
}