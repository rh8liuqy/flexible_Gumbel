rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(coda)
library(MASS)
library(tidyverse)
library(readxl)
library(gridExtra)
library(latex2exp)
library(parallel)
library(HDInterval)
color_scheme_set("brightblue")

## data import
df1 <- readxl::read_xlsx("./crime.xlsx")

y <- df1$`murder rate`
X <- as.matrix(df1[,c("college","poverty","metropolitan")])
N <- nrow(X)
P <- ncol(X)

## data list for stan

dat <- list(N = N,
            P = P,
            X = X,
            y = y)

## FG reg

stan_file <- "../FG_MLR.stan"
stan_mod <- cmdstan_model(stan_file)

## MCMC
stan_fit <- stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  iter_warmup = 10000,
  save_warmup = FALSE,
  iter_sampling = 100000,
)

## parameter estimation 
par.est <- stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
print(par.est)
post_df <- stan_fit$draws(c("alpha",paste0("beta[",1:P,"]")), format = "df")
hdi(post_df$alpha)
hdi(post_df$`beta[1]`)
hdi(post_df$`beta[2]`)
hdi(post_df$`beta[3]`)

#mcmc_trace(post_df)

## normal reg

mean_reg <- lm(df1$`murder rate` ~ df1$college + df1$poverty + df1$metropolitan)
summary(mean_reg)
confint(mean_reg,level = 0.95)

## define negative loglikelihood of FG model
library(stats4)
source("../density.R")

designX <- cbind(1,X)

FG_neglog <- function(w1,sigma1,sigma2,alpha,beta1,beta2,beta3) {
  # w1 <- pars[1]
  # sigma1 <- pars[2]
  # sigma2 <- pars[3]
  beta <- as.vector(c(alpha,beta1,beta2,beta3))
  eta <- as.numeric(designX %*% beta)
  p1 <- w1*rgumbel_pdf(y,eta,sigma1,log=FALSE)
  p2 <- (1-w1)*lgumbel_pdf(y,eta,sigma2,log=FALSE)
  den_value <- p1 + p2
  output <- -sum(log(den_value))
  return(output)
}

FG_mle <- mle(FG_neglog,
              start = list(w1 = 0.9,
                           sigma1 = 2,
                           sigma2 = 52,
                           alpha = -24,
                           beta1 = 0.5,
                           beta2 = 1,
                           beta3 = 0.05),
              lower = c(1e-6,1e-6,1e-6,-Inf,-Inf,-Inf,-Inf),
              upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf))
summary_FG <- summary(FG_mle)
tab_FG <- data.frame(summary_FG@coef)
tab_FG$lower <- tab_FG$Estimate - qnorm(0.975)*tab_FG$Std..Error
tab_FG$upper <- tab_FG$Estimate + qnorm(0.975)*tab_FG$Std..Error
sapply(tab_FG[4:7,], function(x){round(x,3)})
