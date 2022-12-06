## Heritability analysis for YAPP QC'd Dataset
## By John McAuley
## Dec 6th 2022


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
grm.auto <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/workfiles.GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/workfiles.GRM.adj.grm.id")

# Read in pedigree &  ID-Recode-Key
idkey <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ID-Recode-Key.txt", header = T, stringsAsFactors = F)[,c(2, 4)]
names(idkey) <- c("ringnr", "id")
ped <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/SNP_pedigree_Helgeland_05122017.txt", header = TRUE, stringsAsFactors = F)
names(ped)[1] <- "ringnr"
ped <- left_join(ped, idkey)
idkey$id <- as.character(idkey$id)
ped$id <- as.character(ped$id)

#~~ Morph & Brood data
Morph <- read.csv("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/AdultMorphology-pre4_SJ.csv", stringsAsFactors = F)
brood <- read.table("C:/Users/johnb/Projects/PhD_Repo/data/Crimap/Crimap_Troubleshoot_July2021/sparrowabel_crimap/sng_snp_removed/brood.txt", header = TRUE)
brood <- brood[c(2,5)]
names(brood) <- c("ringnr","brood_id")
brood$ringnr <- as.character(brood$ringnr)
#brood$id <- as.character(brood$id)
basedata <- subset(Morph, select = c(ringnr, sex, hatchisland, hatchyear)) %>% unique
basedata$sex <- NULL

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

#  Set character variables to factors
sparrow[sapply(sparrow, is.character)] <- lapply(sparrow[sapply(sparrow, is.character)], as.factor)



#_________________________________________________________________________________________________________
#  Run ANIMAL MODELS
#_________________________________________________________________________________________________________

