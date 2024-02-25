rm(list=ls())
library(rstan)
library(R2jags)
library(coda)
options(mc.cores = parallel::detectCores())

df1 <- read.csv("data_toy.csv")

y <- df1$change
N <- nrow(df1)

dat <- list(N = N,
            y = y)

##rstan
fit <- stan(file = 'model_fitting.stan', 
            data = dat,
            iter = 100000,
            seed = 100,
            chains = 15)

summary(fit)

# AIC of FG ---------------------------------------------------------------

sum_fit <- summary(fit)

sum_fit$summary[,6][1:4]

loglik_value <- sum(mixture_gumbel_pdf(x = df1$change,
                                       w1 = sum_fit$summary[,6][1:4][1],
                                       loc = sum_fit$summary[,6][1:4][2],
                                       scale1 = sum_fit$summary[,6][1:4][3],
                                       scale2 = sum_fit$summary[,6][1:4][4],
                                       log=TRUE))
AIC_value <- -2*loglik_value + 2*(4) ##2506.299

# BIC of FG model ---------------------------------------------------------

BIC_value <- -2*loglik_value + 4*log(length(df1$change)) ##2521.909

#rstan::traceplot(fit)
saveRDS(fit,"stanoutput.rds")

##JAGS
jags_dat <- list(N=N,
                 y=y,
                 ones=rep(1,N))

jagsfit <- jags.parallel(data = jags_dat,
                         parameters.to.save = c("w1","theta","sigma1","sigma2"),
                         model.file = "Gumbel_Mixture_JAGS.BUG",
                         n.chains = 1,
                         n.iter = 10000,
                         n.thin = 1,
                         jags.seed = 100)
print(jagsfit)
mcmcobj <- mcmc(jagsfit$BUGSoutput$sims.matrix[,2:5],
                start=jagsfit$BUGSoutput$n.burnin+1,
                end=jagsfit$BUGSoutput$n.iter,
                thin=jagsfit$BUGSoutput$n.thin)
plot(mcmcobj)