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
source("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ASReml4.EstEffects.R")
load("results/COQCed_models.RData")

wald.asreml(model.basic)
summary.asreml(model.basic, coef = T)$coef.fixed # Fixed effects
summary.asreml(model.basic, coef = T)$varcomp    # random effects
asreml4pin(model.basic)

create_model_table <- function(model, respons.var, mod.struc){

  x <- as.data.frame(wald.asreml(model))
  x$model_term <- row.names(x)
  y <- as.data.frame(summary.asreml(model, coef = T)$coef.fixed) # Fixed effects
  y$model_term <- row.names(y)
  x <- full_join(x,y)
  row.names(x) <- x$model_term
  x$model_term <- NULL
  write.csv(x, file = paste0("results/Model_Tables/",respons.var,"/",mod.struc,"_Fixed.csv"), row.names = TRUE)

  z <- summary.asreml(model, coef = T)$varcomp    # random effects
  z$Effect <- row.names(z)
  p <- asreml4pin(model) # 
  p <- left_join(p,z)
  row.names(p) <- p$Effect
  p$Effect <- NULL
  write.csv(p, file = paste0("results/Model_Tables/",respons.var,"/",mod.struc,"_Random.csv"), row.names = TRUE)
}

create_model_table(model.basic,"Shuff","Basic_AllSpar")
create_model_table(model.basic.F,"Shuff","Basic_Female")
create_model_table(model.basic.M,"Shuff","Basic_Male")
create_model_table(model.basic.F.Mat,"Shuff","Maternal_Female")
create_model_table(model.basic.Mat,"Shuff","Maternal_AllSpar")
create_model_table(model.basic.M.Mat,"Shuff","Maternal_Male")
create_model_table(model.basic.F.Pat,"Shuff","Paternal_Female")
create_model_table(model.basic.Pat,"Shuff","Paternal_AllSpar")
create_model_table(model.basic.M.Pat,"Shuff","Paternal_Male")
create_model_table(model.basic.all,"Shuff","FullModel_AllSpar")
create_model_table(model.basic.F.all,"Shuff","FullModel_Female")
create_model_table(model.basic.M.all,"Shuff","FullModel_Male")