rm(list=ls())
source("../density.R")
source("../EM_v1.0.R")
source("../Bayes_MH_v4.0.r")
library(mixtools)
library(parallel)

## right skewed gumbel mixture
dgumbelmixture_right_skewed <- function(x,w1,loc,scale1,scale2,log=FALSE){
  if (log) {
    if ((w1<0) | (w1>1)) {
      return(-Inf)
    }
    else {
      p1 <- log(w1)+rgumbel_pdf(x,loc,scale1,log=TRUE)
      p2 <- log(1-w1)+rgumbel_pdf(x,loc,scale2,log=TRUE)
      output <- log_sum_exp(p1,p2)
      return(output)
    }
  }
  else {
    if ((w1<0) | (w1>1)) {
      return(0)
    }
    else {
      p1 <- w1*rgumbel_pdf(x,loc,scale1)
      p2 <- (1-w1)*rgumbel_pdf(x,loc,scale2)
      output <- p1+p2
      return(output)
    }
  }
}

## right skewed gumbel generation
rgumbelmixture_right_skewed <- function(n,w1,loc,scale1,scale2){
  n <- n
  Z <- rbinom(n=n,size=1,prob=w1)
  X1 <- rgumbel_RNG(n=n,loc=loc,scale=scale1)
  X2 <- rgumbel_RNG(n=n,loc=loc,scale=scale2)
  Y <- Z*X1 + (1-Z)*X2
  return(Y)
}

simu <- function(seed){
  set.seed(seed)
  y <- rgumbelmixture_right_skewed(200,0.5,0,2,6)
  #sample for KL distance calculation
  KLy <- rgumbelmixture_right_skewed(50000,0.5,0,2,6)
  #ECM HGMD
  est_gumbel <- em(x = y,w1 = 0.5,loc = 0,scale1 = 1, scale2 = 1)
  p <- dgumbelmixture_right_skewed(KLy,w1=0.5,loc=0,scale1=2,scale2=6,log=TRUE)
  q <- mixture_gumbel_pdf(KLy,
                          est_gumbel$w1,
                          est_gumbel$loc,
                          est_gumbel$scale1,
                          est_gumbel$scale2,
                          log=TRUE)
  KLdistance_HGMD_ECM <- mean(p-q)
  #Bayes HGMD
  output <- MH(data=y,
               w1_initial=0.5,
               loc_initial=0,
               scale1_initial=1,
               scale2_initial=1,
               loc_tune=1,
               scale1_tune=1,
               scale2_tune=2)
  est_gumbel <- data.frame(t(summary(output)[[2]][,3]))
  q <- mixture_gumbel_pdf(KLy,
                          est_gumbel$w1,
                          est_gumbel$loc,
                          est_gumbel$scale1,
                          est_gumbel$scale2,
                          log=TRUE)
  KLdistance_HGMD_Bayes <- mean(p-q)
  
  #EM normal mixture
  est_normal <- normalmixEM(y, lambda = .5, mu = c(0.3,4), sigma = c(2,8))
  q <- log(est_normal$lambda[1]*
             dnorm(KLy,
                   mean=est_normal$mu[1],
                   sd=est_normal$sigma[1])+
             est_normal$lambda[2]*
             dnorm(KLy,
                   mean=est_normal$mu[2],
                   sd=est_normal$sigma[2]))
  KLdistance_normal <- mean(p-q)
  df_out <- data.frame(KL_gumbel_ECM=KLdistance_HGMD_ECM,
                       KL_gumbel_Bayes=KLdistance_HGMD_Bayes,
                       KL_normal=KLdistance_normal)
  return(df_out)
}

output <- mclapply(1:1000,
                   simu,
                   mc.cores = 15)

#df_out <- bind_rows(output)

## filter out failed estimation
index <- sapply(output, class) == "data.frame"
output <- output[index]
df_out <- bind_rows(output)
write.csv(df_out,file="right_gumbel_simu.csv",row.names=FALSE)

p1 <- df_out %>%
  gather() %>%
  ggplot(aes(x = key, y = value)) +
  geom_boxplot() +
  ylab("Kullback - Leibler divergence from mixture to right_gumbel") +
  ggtitle("right_gumbel")

p1