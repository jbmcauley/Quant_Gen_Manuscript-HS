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
grm.auto <- read.table("C:/Users/johnb/Projects/Hert_Manuscript-git/data/GWAS/Recreate_Sparrowgen/gcta/gcta_data1_adj_1.94.grm.gz")
ids.auto <- read.table("C:/Users/johnb/Projects/Hert_Manuscript-git/data/GWAS/Recreate_Sparrowgen/gcta/gcta_data1_adj_1.94.grm.id")

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

#Make your GRM
grminv <- makeGRM(grm.auto, ids.auto, id.vector = sparrow$parent) # vector of parents from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

# Basic Structure with Total_coverage as fixed effect
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2)
model.basic <- asreml(fixed = yapp_CO_count_QCed ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual = ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual = ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_QCed, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual = ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

asreml4pin(model.basic)
asreml4pin(model.basic.M)
asreml4pin(model.basic.F)

# Expansion of previous structure with: Maternal Random Effects

#Make your GRM
grm.vec <- unique(c(sparrow$parent, sparrow$Mother))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

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

grm.vec <- unique(c(sparrow$parent, sparrow$Father))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

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

rm(model.basic)
rm(model.basic.F)
rm(model.basic.M)
rm(model.basic.Mat)
rm(model.basic.F.Mat)
rm(model.basic.M.Mat)
rm(model.basic.all)
rm(model.basic.F.all)
rm(model.basic.M.all)
rm(model.basic.Pat)
rm(model.basic.F.Pat)
rm(model.basic.M.Pat)
#____________________________________________________________________________________________
#Intra, Intra-gene, Macro-chrs
#____________________________________________________________________________________________
#Intra Shuff Gene

grm.vec <- unique(c(sparrow$parent))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed)
model.basic <- asreml(fixed = intra_shuff_gene ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Maternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Mother))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed)
model.basic.Mat <- asreml(fixed = intra_shuff_gene ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F.Mat <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M.Mat <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")



# Expansion of previous structure with: Maternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland, yapp_CO_count_QCed)
model.basic.all <- asreml(fixed = intra_shuff_gene ~ 1 + SexF + Total_Coverage + Total_Coverage2 + par_age + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F.all <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + par_age + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M.all <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + par_age + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) + par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

# Expansion of previous structure with: Paternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Father))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Father, yapp_CO_count_QCed)
model.basic.Pat <- asreml(fixed = intra_shuff_gene ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Father, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M.Pat <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff_gene, SexF, Total_Coverage, Total_Coverage2, Father, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F.Pat <- asreml(fixed = intra_shuff_gene ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

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
     ,model.basic.M.Pat, file = "results/intra_shuff_gene_models.RData")

rm(model.basic)
rm(model.basic.F)
rm(model.basic.M)
rm(model.basic.Mat)
rm(model.basic.F.Mat)
rm(model.basic.M.Mat)
rm(model.basic.all)
rm(model.basic.F.all)
rm(model.basic.M.all)
rm(model.basic.Pat)
rm(model.basic.F.Pat)
rm(model.basic.M.Pat)

#Intra_Shuffling

grm.vec <- unique(c(sparrow$parent))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed)
model.basic <- asreml(fixed = intra_shuff ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Maternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Mother))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed)
model.basic.Mat <- asreml(fixed = intra_shuff ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed) %>% filter(SexF == "Female")
model.basic.F.Mat <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, Mother, yapp_CO_count_QCed) %>% filter(SexF == "Male")
model.basic.M.Mat <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")



# Expansion of previous structure with: Maternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Mother, par_age, off_hatchisland, par_hatchisland)
model.basic.all <- asreml(fixed = intra_shuff ~ 1 + SexF + Total_Coverage + Total_Coverage2 +par_age + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Female")
model.basic.F.all <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + par_age + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Male")
model.basic.M.all <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + par_age + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) + par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

# Expansion of previous structure with: Paternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Father))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Father)
model.basic.Pat <- asreml(fixed = intra_shuff ~ 1 + SexF + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Father) %>% filter(SexF == "Male")
model.basic.M.Pat <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, intra_shuff, SexF, Total_Coverage, Total_Coverage2, yapp_CO_count_QCed, Father) %>% filter(SexF == "Female")
model.basic.F.Pat <- asreml(fixed = intra_shuff ~ 1 + Total_Coverage + Total_Coverage2 + yapp_CO_count_QCed,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

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
     ,model.basic.M.Pat, file = "results/intra_shuff_models.RData")

rm(model.basic)
rm(model.basic.F)
rm(model.basic.M)
rm(model.basic.Mat)
rm(model.basic.F.Mat)
rm(model.basic.M.Mat)
rm(model.basic.all)
rm(model.basic.F.all)
rm(model.basic.M.all)
rm(model.basic.Pat)
rm(model.basic.F.Pat)
rm(model.basic.M.Pat)

# Basic Structure with Total_coverage as fixed effect

grm.vec <- unique(c(sparrow$parent))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2)
model.basic <- asreml(fixed = yapp_CO_count_macro ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + ide(parent),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Male")
model.basic.M <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2) %>% filter(SexF == "Female")
model.basic.F <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                        random = ~ vm(parent, grminv) + ide(parent),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Maternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Mother))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother)
model.basic.Mat <- asreml(fixed = yapp_CO_count_macro ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother) %>% filter(SexF == "Female")
model.basic.F.Mat <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother) %>% filter(SexF == "Male")
model.basic.M.Mat <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")



# Expansion of previous structure with: Maternal Random Effects
sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland)
model.basic.all <- asreml(fixed = yapp_CO_count_macro ~ 1 + SexF + Total_Coverage + Total_Coverage2 +par_age,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Female")
model.basic.F.all <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2 + par_age,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother)+ par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Mother, par_age, off_hatchisland, par_hatchisland) %>% filter(SexF == "Male")
model.basic.M.all <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2 + par_age,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Mother, grminv) + ide(Mother) + par_hatchisland + off_hatchisland,
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")




# Expansion of previous structure with: Paternal Random Effects

grm.vec <- unique(c(sparrow$parent, sparrow$Father))
grminv <- makeGRM(grm.auto, ids.auto, id.vector = grm.vec) # vector of IDs from the dataset that you use for the asreml model
# You MUST specify this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Father)
model.basic.Pat <- asreml(fixed = yapp_CO_count_macro ~ 1 + SexF + Total_Coverage + Total_Coverage2,
                          random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Father) %>% filter(SexF == "Male")
model.basic.M.Pat <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                            random = ~ vm(parent, grminv) + ide(parent) + vm(Father, grminv) + ide(Father),
                            residual= ~ idv(units),
                            na.action = na.method(x = "omit", y = "omit"),
                            data = sparrow.sub, workspace = "8gb")

sparrow.sub <- sparrow %>% dplyr::select(parent, yapp_CO_count_macro, SexF, Total_Coverage, Total_Coverage2, Father) %>% filter(SexF == "Female")
model.basic.F.Pat <- asreml(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
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
     ,model.basic.M.Pat, file = "results/macro_models.RData")
rm(model.basic)
rm(model.basic.F)
rm(model.basic.M)
rm(model.basic.Mat)
rm(model.basic.F.Mat)
rm(model.basic.M.Mat)
rm(model.basic.all)
rm(model.basic.F.all)
rm(model.basic.M.all)
rm(model.basic.Pat)
rm(model.basic.F.Pat)
rm(model.basic.M.Pat)


#____________________________________________________________________________________________________
# Write a function for running a model to reduce the clutter above (WIP)
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


#____________________________________________________________________________________________________
# Genetic Correlations (Bivariate Model)
#____________________________________________________________________________________________________

sparrow$MaleACC  <- sparrow$yapp_CO_count_QCed
sparrow$FemaleACC <- sparrow$yapp_CO_count_QCed

# Then remove the values you don't want. One solution is to index by the
# differentiating factor (in this case sex) and change the unwanted values to
# NAs. For example:

sparrow$MaleACC[which(sparrow$SexF == "Female")] <- NA   # if sex == 0 (female) then make NA
sparrow$FemaleACC[which(sparrow$SexF == "Male")] <- NA   # if sex == 1 ( male ) then make NA

#~~ We then have to use the corgh and idh structure...

biv.model <- asreml(fixed  = cbind(MaleACC,FemaleACC) ~ trait:Total_Coverage + trait:Total_Coverage2,   # remove sex!
                 random = ~ corgh(trait):vm(parent, grminv) + idh(trait):ide(parent),
                 residual   = ~ units:us(trait, init = NA),
                 data = sparrow, 
                 maxit = 20,
                 workspace = "8gb")

summary(biv.model)
model5$vparameters

vpredict(biv.model, ra~V1)

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
