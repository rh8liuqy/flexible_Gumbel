rm(list=ls())
source("../density.R")
library(HDInterval)
library(tidyverse)
library(kableExtra)
library(latex2exp)
ECM_est <- read_rds("../table2/ECM_est.rds")
ECM_SE <- read_rds("../table2/ECM_SE.rds")
Bayes_est <- read_rds("../table2/stanoutput.rds")
Bayes_summary <- summary(Bayes_est)
point.est.b <- Bayes_summary$summary[1:4,6] ##posterior median

df1 <- read.csv("../table2/data_toy.csv")
density_est <- density(df1$change,bw="SJ")
df_density <- data.frame(X = density_est$x,
                         value = density_est$y,
                         type = "KDE-SJ")
df_ECM <- data.frame(X = seq(min(df1$change),
                             max(df1$change),
                             length.out = 1000))
df_ECM$value <- mixture_gumbel_pdf(x=df_ECM$X,
                                   w1=ECM_est$w1,
                                   loc=ECM_est$loc,
                                   scale1=ECM_est$scale1,
                                   scale2=ECM_est$scale2)
df_ECM$type <- "ECM"

df_Bayes <- data.frame(X = seq(min(df1$change),
                               max(df1$change),
                               length.out = 1000))
df_Bayes$value <- mixture_gumbel_pdf(x=df_Bayes$X,
                                     w1=point.est.b[1],
                                     loc=point.est.b[2],
                                     scale1=point.est.b[3],
                                     scale2=point.est.b[4])
df_Bayes$type <- "Bayes"

## Normal mixture
est_normal <- readRDS("../table2/normal_mixture_EM.rds")
df_normal <- data.frame(X = seq(min(df1$change),
                                max(df1$change),
                                length.out = 1000))
df_normal$value <- est_normal$lambda[1]*
  dnorm(df_normal$X,
        mean=est_normal$mu[1],
        sd=est_normal$sigma[1])+
  est_normal$lambda[2]*
  dnorm(df_normal$X,
        mean=est_normal$mu[2],
        sd=est_normal$sigma[2])
df_normal$type <- "normal mixture(EM)"

df_plot <- rbind(df_density,df_ECM,df_Bayes,df_normal)
df_plot$type <- factor(df_plot$type,levels = c("KDE-SJ","ECM","Bayes","normal mixture(EM)"))

p1 <- df_plot %>% ggplot(aes(x = X, y = value, linetype = type, color = type)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("density") +
  xlab("daily maximum water elevation change (in feet)") +
  ggtitle("Three difference density estimation for same water elevation change data.")
p1
ggsave("densityplot-1.pdf",p1,height = 5,width = 7)
