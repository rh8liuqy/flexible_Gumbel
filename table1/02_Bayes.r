rm(list=ls())
library(parallel)
source("../density.R")
source("../Bayes_MH_v4.0.r")
simu <- function(seed,n=100) {
  set.seed(seed)
  w1 <- 0.5
  loc <- 0
  scale1 <- 1
  scale2 <- 5
  Z <- rbinom(n=n,size=1,prob=w1)
  X1 <- rgumbel_RNG(n=n,loc,scale1)
  X2 <- 2*loc-rgumbel_RNG(n=n,loc,scale2)
  Y <- Z*X1 + (1-Z)*X2
  output <- MH(data=Y,
               w1_initial=w1,
               loc_initial=loc,
               scale1_initial=scale1,
               scale2_initial=scale2,
               loc_tune=1,
               scale1_tune=0.85,
               scale2_tune=2)
  estimation <- data.frame(t(summary(output)[[2]][,3]))
  SE <- summary(output)[[1]][,2]
  df_out <- data.frame(w1 = estimation$w1,
                       w1_se = SE[1],
                       loc = estimation$loc,
                       loc_se = SE[2],
                       scale1 = estimation$scale1,
                       scale1_se = SE[3],
                       scale2 = estimation$scale2,
                       scale2_se = SE[4])
  
  return(df_out)
}

#sample size 100
output <- mclapply(1:1000,
                   function(x){simu(seed=x,n=100)},
                   mc.cores=15)
df_simu <- bind_rows(output)
write.csv(df_simu,"Bayes_n100.csv",row.names = FALSE)

#sample size 200
output <- mclapply(1:1000,
                   function(x){simu(seed=x,n=200)},
                   mc.cores=15)
df_simu <- bind_rows(output)
write.csv(df_simu,"Bayes_n200.csv",row.names = FALSE)

# revision - add a simulation with sample size 50 -------------------------

simu <- function(seed,n=100) {
  set.seed(seed)
  w1 <- 0.4
  loc <- 1.0
  scale1 <- 1.0
  scale2 <- 1.0
  Z <- rbinom(n=n,size=1,prob=w1)
  X1 <- rgumbel_RNG(n=n,loc,scale1)
  X2 <- 2*loc-rgumbel_RNG(n=n,loc,scale2)
  Y <- Z*X1 + (1-Z)*X2
  output <- MH(data=Y,
               w1_initial=w1,
               loc_initial=loc,
               scale1_initial=scale1,
               scale2_initial=scale2,
               loc_tune=0.2,
               scale1_tune=0.4,
               scale2_tune=0.3)
  estimation <- data.frame(t(summary(output)[[2]][,3]))
  SE <- summary(output)[[1]][,2]
  df_out <- data.frame(w1 = estimation$w1,
                       w1_se = SE[1],
                       loc = estimation$loc,
                       loc_se = SE[2],
                       scale1 = estimation$scale1,
                       scale1_se = SE[3],
                       scale2 = estimation$scale2,
                       scale2_se = SE[4])
  
  return(df_out)
}

output <- mclapply(1:1000,
                   function(x){simu(seed=x,n=50)},
                   mc.cores=8)
df_simu <- bind_rows(output)
write.csv(df_simu,"Bayes_n50.csv",row.names = FALSE)
