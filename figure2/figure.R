rm(list=ls())
library(tidyverse)
df1 <- read.csv("laplace_simu.csv")
colnames(df1) <- c("FG ECM",
                   "FG Bayes",
                   "Normal Mixture Distribution EM")
df1g <- df1 %>%
  gather()
df1g$type <- "E2"
df2 <- read.csv("right_gumbel_simu.csv")
colnames(df2) <- c("FG ECM",
                   "FG Bayes",
                   "Normal Mixture Distribution EM")
df2g <- df2 %>%
  gather()
df2g$type <- "E3"
df3 <- read.csv("Tdist_simu.csv")
colnames(df3) <- c("FG ECM",
                   "FG Bayes",
                   "Normal Mixture Distribution EM")
df3g <- df3 %>%
  gather()
df3g$type <- "E4"
df_plot <- rbind(df1g,df2g,df3g)
df_plot$key <- factor(df_plot$key,
                      levels = c("FG ECM",
                                 "FG Bayes",
                                 "Normal Mixture Distribution EM"))

p1 <- df_plot %>%
  ggplot(aes(x=type,y=value,color=key)) +
  geom_boxplot() +
  ylab("Kullback-Leibler divergence") +
  xlab(NULL) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Kullback-Leibler divergence for experiments (E2) - (E4)")
ggsave("KLdistance-1.pdf",p1,height = 4.8, width = 6)
