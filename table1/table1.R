library(tidyverse)
library(kableExtra)
df_em_50 <- read.csv("EM_n50.csv")
df_em_100 <- read.csv("EM_n100.csv")
df_em_200 <- read.csv("EM_n200.csv")
df_bayes_50 <- read.csv("Bayes_n50.csv")
df_bayes_100 <- read.csv("Bayes_n100.csv")
df_bayes_200 <- read.csv("Bayes_n200.csv")

df_freq <- data.frame(sample_size = rep(c("$n=50$","$n=100$","$n=200$"),each=4),
                      parameter = rep(c("$\\theta$",
                                        "$\\sigma_1$",
                                        "$\\sigma_2$",
                                        "$w_1$"),3),
                      MLE = c(mean(df_em_50$loc),
                              mean(df_em_50$scale1),
                              mean(df_em_50$scale2),
                              mean(df_em_50$w1),
                              mean(df_em_100$loc),
                              mean(df_em_100$scale1),
                              mean(df_em_100$scale2),
                              mean(df_em_100$w1),
                              mean(df_em_200$loc),
                              mean(df_em_200$scale1),
                              mean(df_em_200$scale2),
                              mean(df_em_200$w1)),
                      hat.s.d. = c(mean(df_em_50$loc_se),
                                   mean(df_em_50$scale1_se),
                                   mean(df_em_50$scale2_se),
                                   mean(df_em_50$w1_se),
                                   mean(df_em_100$loc_se),
                                   mean(df_em_100$scale1_se),
                                   mean(df_em_100$scale2_se),
                                   mean(df_em_100$w1_se),
                                   mean(df_em_200$loc_se),
                                   mean(df_em_200$scale1_se),
                                   mean(df_em_200$scale2_se),
                                   mean(df_em_200$w1_se)),
                      mc.se = c(sqrt(var(df_em_50$loc_se)/1000),
                                sqrt(var(df_em_50$scale1_se)/1000),
                                sqrt(var(df_em_50$scale2_se)/1000),
                                sqrt(var(df_em_50$w1_se)/1000),
                                sqrt(var(df_em_100$loc_se)/1000),
                                sqrt(var(df_em_100$scale1_se)/1000),
                                sqrt(var(df_em_100$scale2_se)/1000),
                                sqrt(var(df_em_100$w1_se)/1000),
                                sqrt(var(df_em_200$loc_se)/1000),
                                sqrt(var(df_em_200$scale1_se)/1000),
                                sqrt(var(df_em_200$scale2_se)/1000),
                                sqrt(var(df_em_200$w1_se)/1000)),
                      s.d.1 = c(sd(df_em_50$loc),
                                sd(df_em_50$scale1),
                                sd(df_em_50$scale2),
                                sd(df_em_50$w1),
                                sd(df_em_100$loc),
                                sd(df_em_100$scale1),
                                sd(df_em_100$scale2),
                                sd(df_em_100$w1),
                                sd(df_em_200$loc),
                                sd(df_em_200$scale1),
                                sd(df_em_200$scale2),
                                sd(df_em_200$w1)),
                      bayes.mean = c(mean(df_bayes_50$loc),
                                     mean(df_bayes_50$scale1),
                                     mean(df_bayes_50$scale2),
                                     mean(df_bayes_50$w1),
                                     mean(df_bayes_100$loc),
                                     mean(df_bayes_100$scale1),
                                     mean(df_bayes_100$scale2),
                                     mean(df_bayes_100$w1),
                                     mean(df_bayes_200$loc),
                                     mean(df_bayes_200$scale1),
                                     mean(df_bayes_200$scale2),
                                     mean(df_bayes_200$w1)),
                      bayes.s.d. = c(mean(df_bayes_50$loc_se),
                                     mean(df_bayes_50$scale1_se),
                                     mean(df_bayes_50$scale2_se),
                                     mean(df_bayes_50$w1_se),
                                     mean(df_bayes_100$loc_se),
                                     mean(df_bayes_100$scale1_se),
                                     mean(df_bayes_100$scale2_se),
                                     mean(df_bayes_100$w1_se),
                                     mean(df_bayes_200$loc_se),
                                     mean(df_bayes_200$scale1_se),
                                     mean(df_bayes_200$scale2_se),
                                     mean(df_bayes_200$w1_se)),
                      bayes.mc.se = c(sqrt(var(df_bayes_50$loc_se)/1000),
                                      sqrt(var(df_bayes_50$scale1_se)/1000),
                                      sqrt(var(df_bayes_50$scale2_se)/1000),
                                      sqrt(var(df_bayes_50$w1_se)/1000),
                                      sqrt(var(df_bayes_100$loc_se)/1000),
                                      sqrt(var(df_bayes_100$scale1_se)/1000),
                                      sqrt(var(df_bayes_100$scale2_se)/1000),
                                      sqrt(var(df_bayes_100$w1_se)/1000),
                                      sqrt(var(df_bayes_200$loc_se)/1000),
                                      sqrt(var(df_bayes_200$scale1_se)/1000),
                                      sqrt(var(df_bayes_200$scale2_se)/1000),
                                      sqrt(var(df_bayes_200$w1_se)/1000)),
                      s.d.2 = c(sd(df_bayes_50$loc),
                                sd(df_bayes_50$scale1),
                                sd(df_bayes_50$scale2),
                                sd(df_bayes_50$w1),
                                sd(df_bayes_100$loc),
                                sd(df_bayes_100$scale1),
                                sd(df_bayes_100$scale2),
                                sd(df_bayes_100$w1),
                                sd(df_bayes_200$loc),
                                sd(df_bayes_200$scale1),
                                sd(df_bayes_200$scale2),
                                sd(df_bayes_200$w1))
                      )
df_freq$mc.se <- df_freq$mc.se
df_freq$bayes.mc.se <- df_freq$bayes.mc.se
df_freq[,3:10] <- round(df_freq[,3:10],4)

df_freq$hat.s.d. <- paste0(sprintf("%.4f",df_freq$hat.s.d.),
                           "(",
                           sprintf("%.4f",df_freq$mc.se),
                           ")")

df_freq$bayes.s.d. <- paste0(sprintf("%.4f",df_freq$bayes.s.d.),
                             "(",
                             sprintf("%.4f",df_freq$bayes.mc.se),
                             ")")
df_freq <- df_freq %>%
  dplyr::select(-c("mc.se","bayes.mc.se"))

colnames(df_freq) <- c("sample size",
                       "parameter",
                       "point.est",
                       "$\\widehat{\\text{s.d.}}$",
                       "s.d.",
                       "point.est",
                       "$\\widehat{\\text{s.d.}}$",
                       "s.d.")

kable(df_freq, 
      booktab = TRUE, 
      escape = FALSE,
      caption = "Frequentist and Bayesian inference for experiment (E1) across 1000 Monte Carlo replicates. point.est stands for mean of 1000 point estimations. $\\widehat{\\text{s.d.}}$ and s.d. stand for mean of the corresponding estimated standard deviations and empirical standard deviations respectively. Numbers in parentheses are 100 × (Monte Carlo standard errors) associated with the averages.",
      label = "E1tab") %>%
  collapse_rows(1) %>%
  row_spec(row = 0, align = "c") %>%
  add_header_above(c(" ",
                     " ",
                     "frequentist" = 3,
                     "Bayesian" = 3)) 

tab <- kable(df_freq, 
             booktab = TRUE, 
             escape = FALSE,
             caption = "Frequentist and Bayesian inference for experiment (E1) across 1000 Monte Carlo replicates. point.est stands for mean of 1000 point estimations. $\\widehat{\\text{s.d.}}$ and s.d. stand for mean of the corresponding estimated standard deviations and empirical standard deviations respectively. Numbers in parentheses are Monte Carlo standard errors associated with the averages.",
             label = "E1tab",
             format = 'latex') %>%
  collapse_rows(1) %>%
  row_spec(row = 0, align = "c") %>%
  add_header_above(c(" ",
                     " ",
                     "frequentist" = 3,
                     "Bayesian" = 3))
tab
writeLines(tab, 'tab.tex')
