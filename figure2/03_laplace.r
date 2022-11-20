rm(list=ls())
source("../density.R")
source("../EM_v1.0.R")
source("../Bayes_MH_v4.0.r")
library(mixtools)
library(parallel)

##laplace generation
rlaplace <- function(n,loc,scale){
  X <- rexp(n,rate=1/scale)
  Y <- rexp(n,rate=1/scale)
  Z <- X - Y + loc
  return(Z)
}

##laplace pdf
dlaplace <- function(x,loc,scale,log=FALSE){
  if (log) {
    if (scale <= 0) {
      return(-Inf)
    }
    else {
      p1 <- -log(2*scale)
      p2 <- -abs(x-loc)/scale
      return(p1+p2)
    }
  }
  
  else {
    if (scale <= 0) {
      return(0)
    }
    else {
      p1 <- 1/(2*scale)
      p2 <- exp(-abs(x-loc)/scale)
      output <- p1*p2
      return(output)
    }
  }
}

simu <- function(seed){
  set.seed(seed)
  y <- rlaplace(200,0,2)
  #sample for KL distance calculation
  KLy <- rlaplace(50000,0,2)
  #ECM HGMD
  est_gumbel <- em(x = y,w1 = 0.5,loc = 0,scale1 = 2, scale2 = 2)
  p <- dlaplace(KLy,loc=0,scale=2,log=TRUE)
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
               scale1_initial=2,
               scale2_initial=2,
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
  est_normal <- normalmixEM(y, lambda = .5, mu = c(-0.1,0.1), sigma = c(2.9,3.1))
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
write.csv(df_out,file="laplace_simu.csv",row.names=FALSE)

p1 <- df_out %>%
  gather() %>%
  ggplot(aes(x = key, y = value)) +
  geom_boxplot() +
  ylab("Kullback - Leibler divergence from mixture to laplace") +
  ggtitle("Laplace(0,2)")

p1
