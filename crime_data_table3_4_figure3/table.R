rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(coda)
library(MASS)
library(tidyverse)
library(readxl)
library(gridExtra)
library(latex2exp)
library(parallel)
library(HDInterval)
library(kableExtra)
color_scheme_set("brightblue")

## data import
df1 <- readxl::read_xlsx("./crime.xlsx")

## normal reg

mean_reg <- lm(df1$`murder rate` ~ df1$college + df1$poverty + df1$metropolitan)
sum_tab <- summary(mean_reg)
confint(mean_reg,level = 0.95)

## AIC 303.1543
AIC(mean_reg)
## BIC 312.8134
BIC(mean_reg)

df <- data.frame(sum_tab$coefficients)
df <- df[1:2]
df <- cbind(df,confint(mean_reg,level = 0.95))
df <- apply(df, MARGIN = 2,function(x){round(x,3)})
df
kbl(df,booktabs = TRUE, format = "latex")
