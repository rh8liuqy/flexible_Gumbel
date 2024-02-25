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

# AIC of FG model ---------------------------------------------------------

loglik_value <- sum(mixture_gumbel_pdf(x = df1$change,
                                       w1 = ECM_est$w1,
                                       loc = ECM_est$loc,
                                       scale1 = ECM_est$scale1,
                                       scale2 = ECM_est$scale2,
                                       log=TRUE))
AIC_value <- -2*loglik_value + 2*(4) ##2506.028

# BIC of FG model ---------------------------------------------------------

BIC_value <- -2*loglik_value + 4*log(length(df1$change)) ##2521.638

# AIC of NM model ---------------------------------------------------------

AIC_value <- -2*est_normal$loglik + 2*(5) ##2499.821
BIC_value <- -2*est_normal$loglik + 5*log(length(df1$change)) ##2519.334
