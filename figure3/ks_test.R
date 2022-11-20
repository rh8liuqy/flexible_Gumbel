rm(list=ls())
source("../density.R")
source("../EM_v1.0.R")
source("../Bayes_MH_v4.0.r")
library(tidyverse)
library(goftest)
set.seed(100)

##Right Skewed Gumbel cdf
rgumbel_cdf <- function(x,loc,scale) {
  if (scale > 0) {
    z <- (x-loc)/scale
    output <- exp(-exp(-z))
    return(output)
  }
  else if (scale <= 0){
    print("scale parameter must be > 0!")
    return(0)
  }
}
##Left skewed Gumbel cdf
lgumbel_cdf <- function(x,loc,scale) {
  if (scale > 0) {
    z <- (x-loc)/scale
    output <- 1 - exp(-exp(z))
    return(output)
  }
  else if (scale <= 0){
    print("scale parameter must be > 0!")
    return(0)
  }
}

##mixture cdf
mixture_cdf <- function(x,w1,loc1,loc2,scale1,scale2){
  p1 <- w1*rgumbel_cdf(x,loc1,scale1)
  p2 <- (1-w1)*lgumbel_cdf(x,loc2,scale2)
  output <- p1 + p2
  return(output)
}

df1 <- read.csv("../table2/data_toy.csv")

ECM_est <- read_rds("../table2/ECM_est.rds")
Bayes_est <- read_rds("../table2/stanoutput.rds")
Bayes_summary <- summary(Bayes_est)
est_normal <- readRDS("../table2/normal_mixture_EM.rds")

point.est <- c(ECM_est$w1,
               ECM_est$loc,
               ECM_est$scale1,
               ECM_est$scale2)

point.est.b <- Bayes_summary$summary[1:4,6] ##posterior median

ECM_cdf <- function(x){
  output <- mixture_cdf(x,
                        w1=ECM_est$w1,
                        loc1=ECM_est$loc,
                        loc2=ECM_est$loc,
                        scale1=ECM_est$scale1,
                        scale2=ECM_est$scale2)
  return(output)
}

Bayes_cdf <- function(x){
  output <- mixture_cdf(x,
                        w1=point.est.b[1],
                        loc1=point.est.b[2],
                        loc2=point.est.b[2],
                        scale1=point.est.b[3],
                        scale2=point.est.b[4])
  return(output)
}

Normal_mixture_cdf <- function(x){
  w1 <- est_normal$lambda[1]
  loc1 <- est_normal$mu[1]
  loc2 <- est_normal$mu[2]
  scale1 <- est_normal$sigma[1]
  scale2 <- est_normal$sigma[2]
  p1 <- w1*pnorm(x,loc1,scale1)
  p2 <- (1-w1)*pnorm(x,loc2,scale2)
  return(p1+p2)
}

gumbel_mixture_simu <- function(w1,loc,scale1,scale2,n){
  Z <- rbinom(n=n,size=1,prob=w1)
  X1 <- rgumbel_RNG(n=n,loc,scale1)
  X2 <- 2*loc-rgumbel_RNG(n=n,loc,scale2)
  Y <- Z*X1 + (1-Z)*X2
  return(Y)
}

normal_mixture_simu <- function(w1,loc1,loc2,scale1,scale2,n){
  Z <- rbinom(n=n,size=1,prob=w1)
  X1 <- rnorm(n=n,loc1,scale1)
  X2 <- rnorm(n=n,loc2,scale2)
  Y <- Z*X1 + (1-Z)*X2
  return(Y)
}

## output to show in footnote
df_output <- data.frame(ks_test = rep(0,3))

##KS test with estimated parameter
ks_em_D <- ks.test(df1$change,"ECM_cdf",
                   alternative = "two.sided")
ks_bayes_D <- ks.test(df1$change,"Bayes_cdf",
                      alternative = "two.sided")
ks_normal_mixture_D <- ks.test(df1$change,
                               "Normal_mixture_cdf",
                               alternative = "two.sided")
ks_out <- matrix(0,nrow = 500000, ncol = 3)
for (i in 1:500000) {
  data_ecm <- gumbel_mixture_simu(w1=ECM_est$w1,
                                 loc=ECM_est$loc,
                                 scale1=ECM_est$scale1,
                                 scale2=ECM_est$scale2,
                                 n=nrow(df1))
  data_bayes <- gumbel_mixture_simu(w1=point.est.b[1],
                                    loc=point.est.b[2],
                                    scale1=point.est.b[3],
                                    scale2=point.est.b[4],
                                    n=nrow(df1))
  data_normal_mix <- normal_mixture_simu(w1=est_normal$lambda[1],
                                         loc1=est_normal$mu[1],
                                         loc2=est_normal$mu[2],
                                         scale1=est_normal$sigma[1],
                                         scale2=est_normal$sigma[2],
                                         n=nrow(df1))
  ks_out[i,1] <- ks.test(data_ecm,"ECM_cdf")$statistic
  ks_out[i,2] <- ks.test(data_bayes,"Bayes_cdf")$statistic
  ks_out[i,3] <- ks.test(data_normal_mix,"Normal_mixture_cdf")$statistic
}

df_output$ks_test[1] <- mean(ks_em_D$statistic >= ks_out[,1])
df_output$ks_test[2] <- mean(ks_bayes_D$statistic >= ks_out[,2])
df_output$ks_test[3] <- mean(ks_normal_mixture_D$statistic >= ks_out[,3])

df_output
saveRDS(df_output,"gof.rds")
