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
grm.auto <- read.table("C:/Users/johnb/Projects/Hert_Manuscript-git/data/GWAS/Recreate_Sparrowgen/gcta_GRM/gcta_data1_adj_1.94.grm.gz")
ids.auto <- read.table("C:/Users/johnb/Projects/Hert_Manuscript-git/data/GWAS/Recreate_Sparrowgen/gcta_GRM/gcta_data1_adj_1.94.grm.id")

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
sparrow$par_age <- sparrow$off_hatchyear - sparrow$par_hatchyear # Only 60% of data has an 'Age'

#ID Key
idkey <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ID-Recode-Key.txt", header = T, stringsAsFactors = F)[,c(2, 4)]
idkey[,2] <- as.character(idkey[,2])
names(idkey) <- c("parent", "id")
sparrow$parent <- as.character(sparrow$parent)
sparrow <- left_join(sparrow, idkey)


#  Set character variables to factors
sparrow[sapply(sparrow, is.character)] <- lapply(sparrow[sapply(sparrow, is.character)], as.factor)
sparrow$id <- as.factor(sparrow$id)


#_________________________________________________________________________________________________________
#  Run ANIMAL MODELS
#_________________________________________________________________________________________________________


#____________________________________________________________________________________________________
# Write a function for running a model to reduce the clutter above (WIP)
#____________________________________________________________________________________________________

RR_traits <- c("yapp_CO_count_no_micro", "intra_shuff_no_micro", "intra_shuff_gene_no_micro",
               "yapp_CO_count_QCed", "intra_shuff_gene", "intra_shuff",
               "yapp_CO_count_macro", "intra_shuff_macro", "intra_shuff_gene_macro")

sparrow.sub <- sparrow %>% 
  dplyr::select(parent, yapp_CO_count_macro, 
                SexF, Total_Coverage, Total_Coverage2, 
                Mother, par_age, off_hatchisland, par_hatchisland)


fixed.effects <- "par_age"
random.effects <- c("Mother", "off_hatchisland", "par_hatchisland")


grm.vec <- unique(sparrow$parent)
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the data-set that you use for the asreml model
attr(grminv, which = "INVERSE") <- TRUE # You MUST specify this is an inverted matrix!


#Loop for basic models
model.list <- list() #results list
for(m in RR_traits){
  
  m <- as.character(m)
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula1 <- as.formula(paste0(m, "~ SexF + Total_Coverage + Total_Coverage2"))
  }else{fixed_formula1<-as.formula(paste0(m, "~ SexF + Total_Coverage+Total_Coverage2 + yapp_CO_count_no_micro"))}
  
  sparrow.sub <- sparrow %>% dplyr::select(parent, paste0(m), SexF, Total_Coverage, Total_Coverage2)
  model.all <- asreml(fixed = fixed_formula1,
                    random = ~ vm(parent, grminv) + ide(parent),
                    data = sparrow.sub,
                    na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula2 <- as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2"))
  }else{fixed_formula2<-as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2 + yapp_CO_count_no_micro"))}
  
  sparrow.sub <- sparrow %>% dplyr::select(parent, paste0(m), SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
  model.female <- asreml(fixed = fixed_formula2,
                     random = ~ vm(parent, grminv) + ide(parent),
                     data = sparrow.sub,
                     na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  sparrow.sub <- sparrow %>% dplyr::select(parent, paste0(m), SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
  model.male <- asreml(fixed = fixed_formula2,
                        random = ~ vm(parent, grminv) + ide(parent),
                        data = sparrow.sub,
                        na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  model.list[[length(model.list) + 1]] <- model.all
  model.list[[length(model.list) + 1]] <- model.female
  model.list[[length(model.list) + 1]] <- model.male
}      


#Loop for maternal models
# model.list <- list() #results list, currently adds to model.list object
for(m in RR_traits){
  
  m <- as.character(m)
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula1 <- as.formula(paste0(m, "~ SexF + Total_Coverage + Total_Coverage2"))
  }else{fixed_formula1<-as.formula(paste0(m, "~ SexF + Total_Coverage+Total_Coverage2 + yapp_CO_count_no_micro"))}
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland)
  
  model.all <- asreml(fixed = fixed_formula1,
                      random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                      data = sparrow.sub,
                      na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula2 <- as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2"))
  }else{fixed_formula2<-as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2 + yapp_CO_count_no_micro"))}
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Female")
  model.female <- asreml(fixed = fixed_formula2,
                         random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                         data = sparrow.sub,
                         na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Male")
  model.male <- asreml(fixed = fixed_formula2,
                       random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                       data = sparrow.sub,
                       na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  model.list[[length(model.list) + 1]] <- model.all
  model.list[[length(model.list) + 1]] <- model.female
  model.list[[length(model.list) + 1]] <- model.male
}


#Loop for all effects models
# model.list <- list() #results list, currently adds to model.list object
for(m in RR_traits){
  
  m <- as.character(m)
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula1 <- as.formula(paste0(m, "~ SexF + Total_Coverage + Total_Coverage2 + par_age"))
  }else{fixed_formula1<-as.formula(paste0(m, "~ SexF + Total_Coverage+Total_Coverage2 + yapp_CO_count_no_micro  + par_age"))}
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland)
  
  model.all <- asreml(fixed = fixed_formula1,
                      random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) + off_hatchisland + par_hatchisland,
                      data = sparrow.sub,
                      na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  if(m == "yapp_CO_count_no_micro" | m == "yapp_CO_count_QCed" | m == "yapp_CO_count_macro"){
    fixed_formula2 <- as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2 + par_age"))
  }else{fixed_formula2<-as.formula(paste0(m, "~ Total_Coverage + Total_Coverage2 + yapp_CO_count_no_micro + par_age"))}
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Female")
  model.female <- asreml(fixed = fixed_formula2,
                         random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) off_hatchisland + par_hatchisland,
                         data = sparrow.sub,
                         na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  sparrow.sub <- sparrow %>% 
    dplyr::select(parent, yapp_CO_count_macro, 
                  SexF, Total_Coverage, Total_Coverage2, 
                  Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Male")
  model.male <- asreml(fixed = fixed_formula2,
                       random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                       data = sparrow.sub,
                       na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  model.list[[length(model.list) + 1]] <- model.all
  model.list[[length(model.list) + 1]] <- model.female
  model.list[[length(model.list) + 1]] <- model.male
}

#____________________________________________________________________________________________________
# Genetic Correlations (Bivariate Model) No micro acc
#____________________________________________________________________________________________________

grm.vec <- unique(c(sparrow$parent))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow$MaleACC  <- sparrow$yapp_CO_count_no_micro
sparrow$FemaleACC <- sparrow$yapp_CO_count_no_micro

# Then remove the values you don't want. One solution is to index by the
# differentiating factor (in this case sex) and change the unwanted values to
# NAs. For example:

sparrow$MaleACC[which(sparrow$SexF == "Female")] <- NA   # if sex == 0 (female) then make NA
sparrow$FemaleACC[which(sparrow$SexF == "Male")] <- NA   # if sex == 1 ( male ) then make NA

#~~ We then have to use the corgh and idh structure...

biv.model.ACC_noMicro <- asreml(fixed  = cbind(MaleACC,FemaleACC) ~ trait:Total_Coverage + trait:Total_Coverage2,   # remove sex!
                 random = ~ corgh(trait):vm(parent, grminv) + idh(trait):ide(parent),
                 residual   = ~ units:us(trait, init = NA),
                 data = sparrow, 
                 maxit = 20,
                 workspace = "8gb")

#No micro genome

sparrow$MaleACC  <- sparrow$intra_shuff_no_micro
sparrow$FemaleACC <- sparrow$intra_shuff_no_micro

# Then remove the values you don't want. One solution is to index by the
# differentiating factor (in this case sex) and change the unwanted values to
# NAs. For example:

sparrow$MaleACC[which(sparrow$SexF == "Female")] <- NA   # if sex == 0 (female) then make NA
sparrow$FemaleACC[which(sparrow$SexF == "Male")] <- NA   # if sex == 1 ( male ) then make NA

#~~ We then have to use the corgh and idh structure...

biv.model.ACC_noMicroIS <- asreml(fixed  = cbind(MaleACC,FemaleACC) ~ trait:Total_Coverage + trait:Total_Coverage2,   # remove sex!
                                random = ~ corgh(trait):vm(parent, grminv) + idh(trait):ide(parent),
                                residual   = ~ units:us(trait, init = NA),
                                data = sparrow, 
                                maxit = 20,
                                workspace = "8gb")



sparrow$MaleACC  <- sparrow$intra_shuff_gene_no_micro
sparrow$FemaleACC <- sparrow$intra_shuff_gene_no_micro

# Then remove the values you don't want. One solution is to index by the
# differentiating factor (in this case sex) and change the unwanted values to
# NAs. For example:

sparrow$MaleACC[which(sparrow$SexF == "Female")] <- NA   # if sex == 0 (female) then make NA
sparrow$FemaleACC[which(sparrow$SexF == "Male")] <- NA   # if sex == 1 ( male ) then make NA

#~~ We then have to use the corgh and idh structure...

biv.model.ACC_noMicroISG <- asreml(fixed  = cbind(MaleACC,FemaleACC) ~ trait:Total_Coverage + trait:Total_Coverage2,   # remove sex!
                                  random = ~ corgh(trait):vm(parent, grminv) + idh(trait):ide(parent),
                                  residual   = ~ units:us(trait, init = NA),
                                  data = sparrow, 
                                  maxit = 20,
                                  workspace = "8gb")


#________________________________________________________________________________
# Fit alternative GRM - Mother (Investigation of Germline Restricted Chr)
#________________________________________________________________________________

#Make your GRM
grminv <- makeGRM(grm.auto, ids.auto, id.vector = c(sparrow$Mother, sparrow$parent)) # vector of IDs from the data-set that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE


sparrow.sub <- sparrow %>% dplyr::select(parent, Mother, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2)
model.basic <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(Mother, grminv) + ide(Mother),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent,Mother, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent,Mother, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

asreml4pin(model.basic)
asreml4pin(model.basic.M)
asreml4pin(model.basic.F)



sparrow.sub <- sparrow %>% dplyr::select(parent, Mother, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2)
model.basic <- asreml(fixed = yapp_CO_count_macro ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(Mother, grminv) + ide(Mother),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent,Mother, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent,Mother, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(Mother, grminv) + ide(Mother),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

asreml4pin(model.basic)
asreml4pin(model.basic.M)
asreml4pin(model.basic.F)

