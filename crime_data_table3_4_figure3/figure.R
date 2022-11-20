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
color_scheme_set("brightblue")

## data import
df1 <- readxl::read_xlsx("./crime.xlsx")
df1$State.initial <- replace(df1$State.initial,is.na(df1$State.initial),"")


library(latticeExtra)
p1 <- splom(~df1[c(6,4,9,3)])
p1 <- p1 + layer(panel.text(x,y,labels = df1$State.initial,cex = 0.8))
pdf("scatter_plot_matrix.pdf",height = 8,width = 8)
p1
dev.off()
