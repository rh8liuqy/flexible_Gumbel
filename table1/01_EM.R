rm(list=ls())
library(parallel)
source("../density.R")
source("../EM_v1.0.R")

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
  estimation <- em(Y,w1,loc,scale1,scale2)
  SE <- sandwichvar_se(Y,
                       estimation$w1,
                       estimation$loc,
                       estimation$scale1,
                       estimation$scale2)
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
write.csv(df_simu,"EM_n100.csv",row.names = FALSE)

#sample size 200
output <- mclapply(1:1000,
                   function(x){simu(seed=x,n=200)},
                   mc.cores=15)
df_simu <- bind_rows(output)
write.csv(df_simu,"EM_n200.csv",row.names = FALSE)

df_simu %>%
  gather() %>%
  group_by(key) %>%
  summarise(mean=mean(value))

sd(df_simu$loc)
sd(df_simu$scale1)
sd(df_simu$scale2)
sd(df_simu$w1)

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
  estimation <- em(Y,w1,loc,scale1,scale2)
  SE <- sandwichvar_se(Y,
                       estimation$w1,
                       estimation$loc,
                       estimation$scale1,
                       estimation$scale2)
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
                   mc.cores=8)
df_simu <- bind_rows(output)
df_simu <- df_simu[complete.cases(df_simu),] # remove few estimation with numerical errors
write.csv(df_simu,"EM_n50.csv",row.names = FALSE)

df_simu %>%
  gather() %>%
  group_by(key) %>%
  summarise(mean=mean(value))
