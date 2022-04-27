library(loo)
library(parallel)

source("RichardsModelFunctions.r")

# priors
bmu1 <- c(5, 2, 1.5, 5)
bsd1 <- c(2, 1, 1, 2)
bmu2 <- c(5, -0.5, 1.5, 5)
bsd2 <- c(2, 1, 1, 2)
bmu3 <- c(5, 2, 1.5, 5)
bsd3 <- c(2, 1, 1, 2)
atau <- 0.001; btau <- 0.001

# starting values
tau <- 1/(0.2^2)
betastart1 <- c(bmu1, tau)
betastart2 <- c(bmu2, tau)
betastart3 <- c(bmu3, tau)

iter <- 10000
reps <- 1000
saveid <- 9001:iter

deltas <- c(seq(5, 1.1, by=-0.2), c(0.99, 0.8, 2/3, 0.49, 0.25, 0.01, 0, -0.25, -0.5, -0.75, -1))

res <- mclapply(deltas, function(d){
  wres <- replicate(reps, {
    beta <- c(5, d, 1.5, 5, tau)
    bmu1 <- c(5, max(c(1, d)), 1.5, 5)
    bmu2 <- c(5, min(c(1, d)), 1.5, 5)
    x <- runif(10, 0, 10) 
    y <- rnorm(length(x), richards(x, beta), 1/sqrt(beta[5]))
    ll1 <- matrix(NA, nrow=iter, ncol=length(x))
    ll2 <- matrix(NA, nrow=iter, ncol=length(x))
    ll3 <- matrix(NA, nrow=iter, ncol=length(x))
    for (k in 1:length(x)){
      rmh1 <- richardsMH1(x[-k], y[-k], bmu1, bsd1, atau, btau, betastart1, iter)
      rmh2 <- richardsMH2(x[-k], y[-k], bmu2, bsd2, atau, btau, betastart2, iter)
      rmh3 <- richardsMH3(x[-k], y[-k], bmu3, bsd3, atau, btau, betastart3, iter)
      ll1[,k] <- apply(log(rmh1$lmat), 1, sum)
      ll2[,k] <- apply(log(rmh2$lmat), 1, sum)
      ll3[,k] <- apply(log(rmh3$lmat), 1, sum)
    }
    logliklist <- list(ll1[saveid,], ll2[saveid,], ll3[saveid,])
    nef1 <- relative_eff(exp(ll1[saveid,]), chain_id=rep(1, length(saveid)))
    nef2 <- relative_eff(exp(ll2[saveid,]), chain_id=rep(1, length(saveid)))
    nef3 <- relative_eff(exp(ll3[saveid,]), chain_id=rep(1, length(saveid)))
    neflist <- list(nef1, nef2, nef3)
    ws <- loo_model_weights(logliklist, method="stacking", r_eff_list = neflist, cores=1)
    return(ws)
  })
  return(wres)
}, mc.cores=4)


#load("stacking_simulation_results.Rda")

pm <- simplify2array(lapply(res, function(m){
  ind <- apply(m, 2, function(x) which(x == max(x)))
  c(mean(ind == 1), mean(ind == 2), mean(ind == 3))
}))

library(ggplot2)
pmd <- data.frame(delta=rep(deltas, 3))
pmd$values <- c(pm[1,], pm[2,], pm[3,])
pmd$Models <- factor(rep(c("Richards (delta > 1)", "Richards (-1 < delta < 1)", "Logistic"), each=length(deltas)))
pmd$section <- c("-1 < delta < 1", "delta > 1")[as.numeric(pmd$delta >= 1)+1]
ggplot(subset(pmd, delta != 0), aes(x=delta, y=values, fill=Models)) +
  geom_area() +
  scale_x_continuous(breaks = seq(-1,5,0.5)) +
  facet_grid(. ~ section, scales="free_x", space="free_x") +
  theme_minimal() +
  scale_fill_grey(labels=expression(f[1](delta==2), f[1], f[2])) +
  ylab("Model Weight") +
  xlab(expression(delta))

