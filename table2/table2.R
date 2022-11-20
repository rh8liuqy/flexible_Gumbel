rm(list=ls())
library(HDInterval)
library(tidyverse)
library(kableExtra)
ECM_est <- read_rds("ECM_est.rds")
ECM_SE <- read_rds("ECM_SE.rds")
Bayes_est <- read_rds("stanoutput.rds")
Bayes_summary <- summary(Bayes_est)
point.est.b <- Bayes_summary$summary[1:4,6] ##posterior median
s.d.b <- Bayes_summary$summary[1:4,3] ##posterior standard deviation
Bayes_CI <- hdi(extract(Bayes_est))
Bayes_CI <- bind_rows(Bayes_CI)[1:4,]
df_ECM_tab <- data.frame(Parameter = c("$w_1$",
                                       "$\\theta$",
                                       "$\\sigma_1$",
                                       "$\\sigma_2$"),
                         point.est = c(ECM_est$w1,
                                       ECM_est$loc,
                                       ECM_est$scale1,
                                       ECM_est$scale2),
                         s.d. = c(ECM_SE$w1_se,
                                  ECM_SE$loc_se,
                                  ECM_SE$scale1_se,
                                  ECM_SE$scale2_se),
                         lower.95 = c(ECM_est$w1 + qnorm(0.025)*ECM_SE$w1_se,
                                      ECM_est$loc + qnorm(0.025)*ECM_SE$loc_se,
                                      ECM_est$scale1 + qnorm(0.025)*ECM_SE$scale1_se,
                                      ECM_est$scale2 + qnorm(0.025)*ECM_SE$scale2_se),
                         upper.95 = c(ECM_est$w1 - qnorm(0.025)*ECM_SE$w1_se,
                                      ECM_est$loc - qnorm(0.025)*ECM_SE$loc_se,
                                      ECM_est$scale1 - qnorm(0.025)*ECM_SE$scale1_se,
                                      ECM_est$scale2 - qnorm(0.025)*ECM_SE$scale2_se),
                         point.est.b = point.est.b,
                         s.d.b = s.d.b,
                         ci.l. = Bayes_CI[,1],
                         ci.u. = Bayes_CI[,2])
df_ECM_tab <- df_ECM_tab[c(2,3,4,1),]
df_ECM_tab[,2:9] <- round(df_ECM_tab[,2:9],4)
colnames(df_ECM_tab) <- c("parameter",
                          "point.est",
                          "$\\widehat{\\text{s.d.}}$",
                          "lower 95",
                          "upper 95",
                          "point.est",
                          "$\\widehat{\\text{s.d.}}$",
                          "lower 95",
                          "upper 95"
)
tab <- kable(df_ECM_tab, 
             booktab = TRUE, 
             escape = FALSE,
             caption = "Frequesntist and Bayesian inferences about daily maximum water elevation change in Lake Murray, South Carolina, United States.",
             label = "E2tab",
             row.names = FALSE,
             format = 'latex') %>%
  add_header_above(c(" " = 1,
                     "frequentist" = 4,
                     "Bayesian" = 4))
tab
writeLines(tab, 'tab.tex')
