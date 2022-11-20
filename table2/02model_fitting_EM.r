rm(list=ls())
library(mixtools)
source("../density.R")
source("../EM_v1.0.R")
source("../Bayes_MH_v4.0.r")

df1 <- read.csv("data_toy.csv")

ECM_est <- em(df1$change,w1=0.5,loc=0,scale1=5,scale2=4)
ECM_SE <- sandwichvar_se(df1$change,
                         ECM_est$w1,
                         ECM_est$loc,
                         ECM_est$scale1,
                         ECM_est$scale2)

est_normal <- normalmixEM(df1$change, lambda = .5, mu = c(-0.1,0.1), sigma = c(2.9,3.1))

saveRDS(ECM_est,"ECM_est.rds")
saveRDS(ECM_SE,"ECM_SE.rds")
saveRDS(est_normal,"normal_mixture_EM.rds")
