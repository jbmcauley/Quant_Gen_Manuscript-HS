## Heritability analysis for YAPP QC'd Dataset
## By John McAuley
## Dec 6th 2022
## R Version 4.2.2


# ___________________________________________________________________________
#  Loading Libraries
# ___________________________________________________________________________
library(asreml)
library(tidyverse)
library(kinship2)
library(ggplot2)
library(reshape2)
library(cowplot)
library(glue)
library(dplyr)
# library(purrr)
# library(furrr)
# library(broom)
source("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ASReml4.EstEffects.R")
source("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/makeGRM.R")

#_________________________________________________________________________________________________________
#  Import Data & Set up Variables
#_________________________________________________________________________________________________________
# Set up GRM
grm.auto <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/workfiles.GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/workfiles.GRM.adj.grm.id")


#~~ Morph & Brood data
Morph <- read.csv("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/AdultMorphology-pre4_SJ.csv", stringsAsFactors = F)
brood <- read.table("C:/Users/johnb/Projects/PhD_Repo/data/Crimap/Crimap_Troubleshoot_July2021/sparrowabel_crimap/sng_snp_removed/brood.txt", header = TRUE)
brood <- brood[c(2,5)]
names(brood) <- c("ringnr","brood_id")
brood$ringnr <- as.character(brood$ringnr)
#brood$id <- as.character(brood$id)
basedata <- subset(Morph, select = c(ringnr, hatchisland, hatchyear)) %>% unique

#Load in Cleaned YAPP Data
load("data/yapp_data_QCed/5_Full_QC_Sparrow_Yapp_Data.RData") #Bring in: linkmap, recchr, recsumm, and rectab
#  linkmap is the yapp linkage map, 
#  recchr is xover count per chr for each meiosis, 
#  recsumm grouped recchr by meiosis, 
#  rectab contains information on the left and right position of xovers.

#  Set up parent and offspring hatch year + Hatch Island
names(basedata) <- c("parent", "par_hatchisland", "par_hatchyear")
sparrow <- left_join(recsumm, basedata, by = "parent")
names(basedata) <- c("offspring", "off_hatchisland", "off_hatchyear")
sparrow <- left_join(sparrow, basedata, by = "offspring")
names(brood) <- c("offspring", "off_brood_id")
sparrow <- left_join(sparrow, brood, by = "offspring")
sparrow$par_age <- sparrow$off_hatchyear - sparrow$par_hatchyear # Only 60% of data has an 'Age' yikes...

#Make your GRM
grminv <- makeGRM(grm.auto, ids.auto, id.vector = sparrow$id) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE
rm(grm.auto)
rm(ids.auto)

#  Set character variables to factors
sparrow[sapply(sparrow, is.character)] <- lapply(sparrow[sapply(sparrow, is.character)], as.factor)



#_________________________________________________________________________________________________________
#  Run ANIMAL MODELS
#_________________________________________________________________________________________________________

#Most basic Models
# sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF)
# model.basic <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF,
#                       random = ~ vm(parent, grminv) + ide(parent),
#                       residual= ~ idv(units),
#                       na.action = na.method(x = "omit", y = "omit"),
#                       data = sparrow.sub, workspace = "8gb")
# 
# sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF) %>% filter(SexF == "Female")
# model.basic.female <- asreml(fixed = yapp_CO_count_QCed ~ 1,
#                       random = ~ vm(parent, grminv) + ide(parent),
#                       residual= ~ idv(units),
#                       na.action = na.method(x = "omit", y = "omit"),
#                       data = sparrow.sub, workspace = "8gb")
# 
# sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF) %>% filter(SexF == "Male")
# model.basic.male <- asreml(fixed = yapp_CO_count_QCed ~ 1,
#                              random = ~ vm(parent, grminv) + ide(parent),
#                              residual= ~ idv(units),
#                              na.action = na.method(x = "omit", y = "omit"),
#                              data = sparrow.sub, workspace = "8gb")
# 



# Expansion of Basic Structure with Total_coverage as fixed effect
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2)
model.basic <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Maternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother)
model.basic.Mat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother) %>% filter(SexF == "Female")
model.basic.F.Mat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother) %>% filter(SexF == "Male")
model.basic.M.Mat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")



# Expansion of previous structure with: Maternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland)
model.basic.all <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2 +par_age,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Female")
model.basic.F.all <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2 + par_age,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Male")
model.basic.M.all <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2 + par_age,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) + par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Paternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Father)
model.basic.Pat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Father) %>% filter(SexF == "Male")
model.basic.M.Pat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2, Father) %>% filter(SexF == "Female")
model.basic.F.Pat <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")


# asreml4pin(model.basic)
# asreml4pin(model.basic.F)
# asreml4pin(model.basic.M)
# asreml4pin(model.basic.Mat)
# asreml4pin(model.basic.F.Mat)
# asreml4pin(model.basic.M.Mat)
# asreml4pin(model.basic.all)
# asreml4pin(model.basic.F.all)
# asreml4pin(model.basic.M.all)
# asreml4pin(model.basic.Pat)
# asreml4pin(model.basic.F.Pat)
# asreml4pin(model.basic.M.Pat)

save(model.basic
,model.basic.F
,model.basic.M
,model.basic.Mat
,model.basic.F.Mat
,model.basic.M.Mat
,model.basic.all
,model.basic.F.all
,model.basic.M.all
,model.basic.Pat
,model.basic.F.Pat
,model.basic.M.Pat, file = "results/COQCed_models.RData")


#____________________________________________________________________________________________________
# Write a function for running a model to reduce the clutter above
#____________________________________________________________________________________________________
var.fixed = c(1, "SexF", "Total_Coverage", "Total_Coverage2")
var.random = c("parent", "Mother")
var.resonse = ("yapp_CO_count_macro")

run_model  <- function(var.response, var.fixed, var.random, dataset){
  asreml(fixed = reformulate(termlabels = var.fixed, response = var.response),
         random = ~,
         residual = ~,
         na.action = na.method(x = "omit", y = "omit"),
         data = dataset,
         workspace = "8gb"
  )
  
}

