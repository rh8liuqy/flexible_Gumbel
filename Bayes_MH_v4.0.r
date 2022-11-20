library(rstan)
library(tidyverse)
library(coda)
library(HDInterval)
library(parallel)

#source("density.R")

MH <- function(data,
               w1_initial=0.5,
               loc_initial=0,
               scale1_initial=1,
               scale2_initial=1,
               loc_tune=1,
               scale1_tune=1,
               scale2_tune=1,
               random_seed=100,
               burn_in=5000,
               iteration=10000,
               thin=1){
  set.seed(random_seed)
  Y <- data
  n <- length(Y)
  w1t <- w1_initial
  loct <- loc_initial
  scale1t <- scale1_initial
  scale2t <- scale2_initial
  zt <- numeric(n)
  pzt <- numeric(n)
  output <- list()
  for (i in 1:iteration) {
    ##Latent Variable Z
    pzt <- (w1t*rgumbel_pdf(Y,loct,scale1t))/
      (w1t*rgumbel_pdf(Y,loct,scale1t)+(1-w1t)*
         lgumbel_pdf(Y,loct,scale2t))
    zt <- sapply(pzt, function(x){rbinom(1,1,x)})
    
    ##w1 prior beta(1,1)
    w1t <- rbeta(1,sum(zt)+1,n-sum(zt)+1)
    
    #loc
    loct_logdensity <- function(x){
      output <- sum(mixture_gumbel_pdf(Y,w1t,x,scale1t,scale2t,log=TRUE))
      return(output)
    }
    loct_can <- rnorm(1,mean=loct,sd=loc_tune)
    r_loc <- exp(loct_logdensity(loct_can)-
                     loct_logdensity(loct))
    r_loc <- min(c(1,r_loc))
    Z_loc <- rbinom(1,1,r_loc)
    loct <- Z_loc*loct_can+(1-Z_loc)*loct
    
    #scale1
    scale1t_logdensity <- function(x){
      if (x <= 0) {
        output <- -Inf
      }
      else {
        output <- sum(mixture_gumbel_pdf(Y,w1t,loct,x,scale2t,log=TRUE)) +
          dinvgamma(x,1,1,log=TRUE)
      }
      return(output)
    }
    scale1t_can <- rnorm(1,scale1t,scale1_tune)
    r_s1 <- exp(scale1t_logdensity(scale1t_can)-
                  scale1t_logdensity(scale1t))
    r_s1 <- min(c(1,r_s1))
    Z_s1 <- rbinom(1,1,r_s1)
    scale1t <- Z_s1*scale1t_can+(1-Z_s1)*scale1t
    
    scale2t_logdensity <- function(x){
      if (x <= 0) {
        output <- -Inf
      }
      else{
      output <- sum(mixture_gumbel_pdf(Y,w1t,loct,scale1t,x,log=TRUE)) +
        dinvgamma(x,1,1,log=TRUE)
      }
      return(output)
    }
    scale2t_can <- rnorm(1,scale2t,scale2_tune)
    r_s2 <- exp(scale2t_logdensity(scale2t_can)-
                  scale2t_logdensity(scale2t))
    r_s2 <- min(c(1,r_s2))
    Z_s2 <- rbinom(1,1,r_s2)
    scale2t <- Z_s2*scale2t_can+(1-Z_s2)*scale2t
    
    output[[i]] <- data.frame(w=w1t,
                              loc=loct,
                              scale1=scale1t,
                              scale2=scale2t)
  }
  df_MCMC <- bind_rows(output)
  df_MCMC <- df_MCMC[(burn_in+1):iteration,]
  sims <- array(NA,
                dim = c(nrow(df_MCMC),1,4), 
                dimnames = list(NULL,
                                NULL,
                                c("w1","loc","scale1","scale2")))
  sims[1:nrow(df_MCMC),1,] <- as.matrix(df_MCMC)
  mcmcobj <- mcmc(sims[1:nrow(df_MCMC),1,],
                  start=burn_in+1,
                  end=iteration,
                  thin=thin)
  return(mcmcobj)
}
