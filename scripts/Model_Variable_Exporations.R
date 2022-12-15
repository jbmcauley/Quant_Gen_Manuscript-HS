## Exploring Models and Creating Tables
## By John McAuley
## Dec 15th 2022
## R Version 4.2.2

library(asreml)
library(tidyverse)
library(kinship2)
library(ggplot2)
library(reshape2)
library(cowplot)
library(glue)
library(dplyr)

load("results/COQCed_models.RData")

wald.asreml(model.basic)
summary.asreml(anmodd3, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd3, coef = T)$varcomp    # random effects
asreml4pin(anmodd3) #

