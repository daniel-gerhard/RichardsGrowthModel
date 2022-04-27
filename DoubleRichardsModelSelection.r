library(loo)

# import model functions and example data
source("DoubleRichardsModelFunctions.r")
grapedata <- read.csv("grapedata.csv")

x <- grapedata$time
y <- log(grapedata$weight)


# truncated normal priors (mean mu and standard deviation sd) for 9 models
bmu1 <- c(4, 2, 0.5, 15, 0.5, 2, 2, 50)
bsd1 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu2 <- c(4, -0.5, 0.5, 15, 0.5, -0.5, 2, 50)
bsd2 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu3 <- c(4, 2, 0.5, 15, 0.5, -0.5, 2, 50)
bsd3 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu4 <- c(4, -0.5, 0.5, 15, 0.5, 2, 2, 50)
bsd4 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu5 <- c(4, 2, 0.5, 15, 0.5, 2, 2, 50)
bsd5 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu6 <- c(4, 2, 0.5, 15, 0.5, 2, 2, 50)
bsd6 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu7 <- c(4, 2, 0.5, 15, 0.5, -0.5, 2, 50)
bsd7 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu8 <- c(4, -0.5, 0.5, 15, 0.5, 2, 2, 50)
bsd8 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
bmu9 <- c(4, 2, 0.5, 15, 0.5, 2, 2, 50)
bsd9 <- c(0.5, 1, 0.5, 2, 0.5, 1, 1, 5)
# inverse gamma priors
atau <- 0.001; btau <- 0.001

iter <- 10000
# 9000 samples warmup
saveid <- 9001:iter

# starting values
betastart1 <- c(bmu1, 0.5)
betastart2 <- c(bmu2, 0.5)
betastart3 <- c(bmu3, 0.5)
betastart4 <- c(bmu4, 0.5)
betastart5 <- c(bmu5, 0.5)
betastart6 <- c(bmu6, 0.5)
betastart7 <- c(bmu7, 0.5)
betastart8 <- c(bmu8, 0.5)
betastart9 <- c(bmu9, 0.5)


# objects to save log-likelihoods
ll1 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll2 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll3 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll4 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll5 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll6 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll7 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll8 <- matrix(NA, nrow=length(saveid), ncol=length(x))
ll9 <- matrix(NA, nrow=length(saveid), ncol=length(x))

# sampling
library(progress)
pb <- progress_bar$new(total = length(x))
for (k in 1:length(x)){
  rmh1 <- richardsMH(x[-k], y[-k], bmu1, bsd1, atau, btau, betastart1, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
  rmh2 <- richardsMH(x[-k], y[-k], bmu2, bsd2, atau, btau, betastart2, iter, ltd1=-1, utd1=1, ltd2=-1, utd2=1, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
  rmh3 <- richardsMH(x[-k], y[-k], bmu3, bsd3, atau, btau, betastart3, iter, ltd1=1, utd1=Inf, ltd2=-1, utd2=1, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
  rmh4 <- richardsMH(x[-k], y[-k], bmu4, bsd4, atau, btau, betastart4, iter, ltd1=-1, utd1=1, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
  rmh5 <- richardsMH(x[-k], y[-k], bmu5, bsd5, atau, btau, betastart5, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, 2, NA, NA, NA, NA, NA, NA))
  rmh6 <- richardsMH(x[-k], y[-k], bmu6, bsd6, atau, btau, betastart6, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, 2, NA, NA))
  rmh7 <- richardsMH(x[-k], y[-k], bmu7, bsd7, atau, btau, betastart7, iter, ltd1=1, utd1=Inf, ltd2=-1, utd2=1, fixed=c(NA, 2, NA, NA, NA, NA, NA, NA))
  rmh8 <- richardsMH(x[-k], y[-k], bmu8, bsd8, atau, btau, betastart8, iter, ltd1=-1, utd1=1, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, 2, NA, NA))  
  rmh9 <- richardsMH(x[-k], y[-k], bmu9, bsd9, atau, btau, betastart9, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, 2, NA, NA, NA, 2, NA, NA))  
  ll1[,k] <- apply(log(rmh1$lmat[saveid,]), 1, sum)
  ll2[,k] <- apply(log(rmh2$lmat[saveid,]), 1, sum)
  ll3[,k] <- apply(log(rmh3$lmat[saveid,]), 1, sum)
  ll4[,k] <- apply(log(rmh4$lmat[saveid,]), 1, sum)
  ll5[,k] <- apply(log(rmh5$lmat[saveid,]), 1, sum)
  ll6[,k] <- apply(log(rmh6$lmat[saveid,]), 1, sum)
  ll7[,k] <- apply(log(rmh7$lmat[saveid,]), 1, sum)
  ll8[,k] <- apply(log(rmh8$lmat[saveid,]), 1, sum)
  ll9[,k] <- apply(log(rmh9$lmat[saveid,]), 1, sum)
  pb$tick()
}
logliklist <- list(ll1, ll2, ll3, ll4, ll5, ll6, ll7, ll8, ll9)
nef1 <- relative_eff(exp(ll1), chain_id=rep(1, length(saveid)))
nef2 <- relative_eff(exp(ll2), chain_id=rep(1, length(saveid)))
nef3 <- relative_eff(exp(ll3), chain_id=rep(1, length(saveid)))
nef4 <- relative_eff(exp(ll4), chain_id=rep(1, length(saveid)))
nef5 <- relative_eff(exp(ll5), chain_id=rep(1, length(saveid)))
nef6 <- relative_eff(exp(ll6), chain_id=rep(1, length(saveid)))
nef7 <- relative_eff(exp(ll7), chain_id=rep(1, length(saveid)))
nef8 <- relative_eff(exp(ll8), chain_id=rep(1, length(saveid)))
nef9 <- relative_eff(exp(ll9), chain_id=rep(1, length(saveid)))
neflist <- list(nef1, nef2, nef3, nef4, nef5, nef6, nef7, nef8, nef9)

# object with model weights
ws <- loo_model_weights(logliklist, method="stacking", r_eff_list = neflist, cores=1)

################################################

# example evaluation for the full data

rmh1 <- richardsMH(x, y, bmu1, bsd1, atau, btau, betastart1, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
rmh2 <- richardsMH(x, y, bmu2, bsd2, atau, btau, betastart2, iter, ltd1=-1, utd1=1, ltd2=-1, utd2=1, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
rmh3 <- richardsMH(x, y, bmu3, bsd3, atau, btau, betastart3, iter, ltd1=1, utd1=Inf, ltd2=-1, utd2=1, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
rmh4 <- richardsMH(x, y, bmu4, bsd4, atau, btau, betastart4, iter, ltd1=-1, utd1=1, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, NA, NA, NA))
rmh5 <- richardsMH(x, y, bmu5, bsd5, atau, btau, betastart5, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, 2, NA, NA, NA, NA, NA, NA))
rmh6 <- richardsMH(x, y, bmu6, bsd6, atau, btau, betastart6, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, 2, NA, NA))
rmh7 <- richardsMH(x, y, bmu7, bsd7, atau, btau, betastart7, iter, ltd1=1, utd1=Inf, ltd2=-1, utd2=1, fixed=c(NA, 2, NA, NA, NA, NA, NA, NA))
rmh8 <- richardsMH(x, y, bmu8, bsd8, atau, btau, betastart8, iter, ltd1=-1, utd1=1, ltd2=1, utd2=Inf, fixed=c(NA, NA, NA, NA, NA, 2, NA, NA))  
rmh9 <- richardsMH(x, y, bmu9, bsd9, atau, btau, betastart9, iter, ltd1=1, utd1=Inf, ltd2=1, utd2=Inf, fixed=c(NA, 2, NA, NA, NA, 2, NA, NA))  


# prediction at days xc
xc <- seq(0, 110, length=250)

# prediction posterior median
pms1 <- rmh1$pm[saveid,]
psamp1 <- sapply(1:nrow(pms1), function(i){
  drichards(xc, pms1[i,])
})
qpr1 <- apply(psamp1, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms2 <- rmh2$pm[saveid,]
psamp2 <- sapply(1:nrow(pms2), function(i){
  drichards(xc, pms2[i,])
})
qpr2 <- apply(psamp2, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms3 <- rmh3$pm[saveid,]
psamp3 <- sapply(1:nrow(pms3), function(i){
  drichards(xc, pms3[i,])
})
qpr3 <- apply(psamp3, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms4 <- rmh4$pm[saveid,]
psamp4 <- sapply(1:nrow(pms4), function(i){
  drichards(xc, pms4[i,])
})
qpr4 <- apply(psamp4, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms5 <- rmh5$pm[saveid,]
psamp5 <- sapply(1:nrow(pms5), function(i){
  drichards(xc, pms5[i,])
})
qpr5 <- apply(psamp5, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms6 <- rmh6$pm[saveid,]
psamp6 <- sapply(1:nrow(pms6), function(i){
  drichards(xc, pms6[i,])
})
qpr6 <- apply(psamp6, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms7 <- rmh7$pm[saveid,]
psamp7 <- sapply(1:nrow(pms7), function(i){
  drichards(xc, pms7[i,])
})
qpr7 <- apply(psamp7, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms8 <- rmh8$pm[saveid,]
psamp8 <- sapply(1:nrow(pms8), function(i){
  drichards(xc, pms8[i,])
})
qpr8 <- apply(psamp8, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pms9 <- rmh9$pm[saveid,]
psamp9 <- sapply(1:nrow(pms9), function(i){
  drichards(xc, pms9[i,])
})
qpr9 <- apply(psamp9, 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))


# weights from model selection procedure
#ws <- c(0.244, 0.000, 0.292, 0.175, 0.000, 0.247, 0.004, 0.039, 0.000) 
#model averaged prediction
qprw <- ws[1]*qpr1 + ws[2]*qpr2 + ws[3]*qpr3 + ws[4]*qpr4 + ws[5]*qpr5 + ws[6]*qpr6 + ws[7]*qpr7 + ws[8]*qpr8 + ws[9]*qpr9 

# prediction for each single model
pdat <- data.frame(xc, 
                   pred=exp(c(qpr1[3,], qpr2[3,], qpr3[3,], qpr4[3,], 
                              qpr5[3,], qpr6[3,], qpr7[3,], qpr8[3,], qpr9[3,])),
                   w=rep(ws, each=length(xc)),
                   model=as.factor(rep(1:9, each=length(xc))))

# model averaged prediction
pwdat <- data.frame(xc, 
                    predm=exp(qprw[3,]),
                    predl=exp(qprw[1,]),
                    predu=exp(qprw[5,]))

library(ggplot2)
ggplot(dat, aes(x=time, y=weight)) +
  geom_point() +
  theme_bw() +
  geom_line(data=pdat, aes(x=xc, y=pred, group=model), colour="grey2") +
  geom_line(data=pwdat, aes(x=xc, y=predm, group=NULL), colour="red2", size=1.05) +
  ylim(c(0, 300)) +
  ylab("Bunch weight [g]") + xlab("Time [days]")



