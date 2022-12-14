## GWAS analysis for YAPP QC'd Dataset
## By John McAuley
## Dec 12th 2022
## R version 3.5.3

# This script is for running the remainder of GWAS models on a server

source("rGLSadj.R")
library(GenABEL)
library(GenABEL.data)
library(RepeatABEL)


#Define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])
print(resultpath)


# Load in phentype data
sparrow <- read.table("sparrow.txt", stringsAsFactors = TRUE, header = TRUE)
sparrow$id <- as.factor(sparrow$id)
sparrow.f <- sparrow[which(sparrow$SexF == "Female"),]
sparrow.m <- sparrow[which(sparrow$SexF == "Male"),]

#Load genabel obj
load("data1.RData")

# gkin setup
load("sparrow.gkin.RData")
sparrow.gkin.sym <- sparrow.gkin
sparrow.gkin.sym[upper.tri(sparrow.gkin.sym)] = t(sparrow.gkin.sym)[upper.tri(sparrow.gkin.sym)]
sparrow.gkin.sym <- sparrow.gkin.sym * 2
rm(sparrow.gkin)


# All sparrows
gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_macro ~ SexF + Total_Coverage + Total_Coverage2, 
                               random = ~1|id,
                               id.name = "id", #YOUR ID MUST BE NAMED "id" otherwise library breaks
                               genabel.data = data1,
                               phenotype.data = sparrow, 
                               corStruc = list(id = list("GRM")),
                               GRM = sparrow.gkin.sym)

gwasACC2 <- rGLSadj(yapp_CO_count_macro ~ SexF + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym)

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)

save(gwasACCres2, file = paste0(resultpath, "AS_Macro-basic_GWAS-lambda-corrected.RData"))
save(gwasACC2, file = paste0(resultpath, "AS_Macro-basic_GWAS.RData"))

rm(gwasACC_prefit2)
rm(gwasACC2)
rm(gwasACCres2)

# Female
gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2, 
                               random = ~1|id,
                               id.name = "id", #YOUR ID MUST BE NAMED "id" otherwise library breaks
                               genabel.data = data1,
                               phenotype.data = sparrow.f, 
                               corStruc = list(id = list("GRM")),
                               GRM = sparrow.gkin.sym)

gwasACC2 <- rGLSadj(yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow.f, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym)

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)

save(gwasACCres2, file = paste0(resultpath, "F_Macro-basic_GWAS-lambda-corrected.RData"))
save(gwasACC2, file = paste0(resultpath, "F_Macro-basic_GWAS.RData"))


# Male
gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2, 
                               random = ~1|id,
                               id.name = "id", #YOUR ID MUST BE NAMED "id" otherwise library breaks
                               genabel.data = data1,
                               phenotype.data = sparrow.m, 
                               corStruc = list(id = list("GRM")),
                               GRM = sparrow.gkin.sym)

gwasACC2 <- rGLSadj(yapp_CO_count_macro ~ 1 + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow.m, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym)

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)

save(gwasACCres2, file = paste0(resultpath, "M_Macro-basic_GWAS-lambda-corrected.RData"))
save(gwasACC2, file = paste0(resultpath, "M_Macro-basic_GWAS.RData"))