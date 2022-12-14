## GWAS analysis for YAPP QC'd Dataset
## By John McAuley
## Dec 12th 2022
## R version 3.5.3

# Since we have repeated measures we are going to use the R package RepeatABEL.
# RepeatABEL also allows for us to account for any underlying population structure
# in our data.
# Notably, this package is no longer actively updated and thus an older version
# of R and its dependencies is typically required in order for it to work. This 
# can cause some annoyance upfront with installation, but once set up RepeatABEL
# is fairly user-friendly. 
# *When analysing large datasets RepeatABEL can be somewhat slow.

# Rtools: You must have a version of Rtools compatible with genabel installation
# versions 4 or greater will likely not work.


#_______________________________________________________
# 0. RepeatABEL installation:
#_______________________________________________________

install.packages("devtools")
library(devtools)

install_github("GenABEL-Project/GenABEL.data")
install_github("GenABEL-Project/GenABEL")
install_github("cran/RepeatABEL")

# Once installed check libraries load without issue:
library(GenABEL)
library(GenABEL.data)
library(RepeatABEL)

# If installation via github fails. One option for installation is through the
# cran archive as follows:

# GenABEL.data
# Download packages from CRAN archive
install.packages("GenABEL.data_1.0.0.tar.gz", repos = NULL)
install.packages("GenABEL_1.8-0.tar.gz", repos = NULL)

install.packages("hglm")
#Download repeatable from archive
install.packages("RepeatABEL_1.1.tar.gz", repos = NULL)

install.packages("dplyr_0.8.0.1.tar.gz", repos = NULL)
#_______________________________________________________
# 1.  Setup 
#_______________________________________________________
# Requirements: Phenotype Data, Genetic data, and GRM
# load rGLSadj.R as it contains functions we will use to extract results
source("C:/Users/johnb/Projects/PhD_Repo/data/YAPP/rGLSadj.R")
library(GenABEL)
library(GenABEL.data)
library(RepeatABEL)
library(dplyr)

# Load in a dataframe with Phenotypic data
sparrow <- read.table("data/yapp_data_QCed/5_Full_Recombination_Phenotypes.txt", header = TRUE)

# Here I load in an idkey b/c my grm has different ids than my phenotype txt file.
idkey <- read.table("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ID-Recode-Key.txt", header = T, stringsAsFactors = F)[,c(2, 4)]
idkey[,2] <- as.character(idkey[,2])
names(idkey) <- c("parent", "id")
sparrow$parent <- as.character(sparrow$parent)
sparrow <- left_join(sparrow, idkey)

# Load in genotype data
load("C:/Users/johnb/Projects/PhD_Repo/data/YAPP/sparrowgen.RData")

# Load GRM
load("C:/Users/johnb/Projects/PhD_Repo/data/YAPP/sparrow.gkin.RData")

# gkin setup
sparrow.gkin.sym <- sparrow.gkin
sparrow.gkin.sym[upper.tri(sparrow.gkin.sym)] = t(sparrow.gkin.sym)[upper.tri(sparrow.gkin.sym)]
sparrow.gkin.sym <- sparrow.gkin.sym * 2
rm(sparrow.gkin)

#Inspect genotype object
qc0snp <- summary.snp.data(sparrowgen@gtdata)
hist(qc0snp$CallRate, breaks = 50)
qc0id <- perid.summary(sparrowgen)
hist(qc0id$CallPP, breaks = 50) 

#Choose Quality Control thresholds for check.marker() function based on histograms
qc1 <- check.marker(sparrowgen, callrate = 0.9, perid.call = 0.8, p.level = 0)
data1 <- sparrowgen[qc1$idok, qc1$snpok]
# save(data1, file = "data/data1.RData")
rm(sparrowgen)
rm(qc1)
rm(qc0id)
rm(qc0snp)

#_______________________________________________________
# 2.  RepeatABEL 
#_______________________________________________________

# Basic Model
gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_QCed ~ SexF + Total_Coverage + Total_Coverage2, 
                               random = ~1|id,
                               id.name = "id", #YOUR ID MUST BE NAMED "id" otherwise library breaks
                               genabel.data = data1,
                               phenotype.data = sparrow, 
                               corStruc = list(id = list("GRM")),
                               GRM = sparrow.gkin.sym)

system.time(gwasACC2 <- rGLSadj(yapp_CO_count_QCed ~ SexF + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)
head(gwasACCres2)
plot_rGLSadj_results(gwasACCres2) +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24,face = "bold"),
        plot.title = element_text(size = 28),
        #plot.margin = margin(20, 50, 20, 50),
        legend.text=element_text(size=24),
        legend.title=element_blank(),
        axis.text.x = element_blank())+
  ggtitle("Both Sexes")+
  scale_colour_grey(start = 0, end = .45)


save(gwasACCres2, file = "results/GWAS/CO_count_QCed/All-sparrow-GWAS-lambda-corrected.RData")
save(gwasACC2, file = "results/GWAS/CO_count_QCed/All-sparrow-GWAS.RData")





# Female Only - Basic Model
sparrow.f <- sparrow %>% filter(SexF == "Female")

system.time(gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2, 
                                           id.name = "id",
                                           genabel.data = data1,
                                           phenotype.data = sparrow.f, 
                                           corStruc = list(id = list("GRM")),
                                           GRM = sparrow.gkin.sym))

system.time(gwasACC2 <- rGLSadj(yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow.f, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)
head(gwasACCres2)
plot_rGLSadj_results(gwasACCres2) +
  theme(axis.text = element_text(size = 24), 
        axis.title = element_text(size = 24,face = "bold"), 
        plot.title = element_text(size = 28),
        #plot.margin = margin(20, 50, 20, 50),
        legend.text=element_text(size=24),
        legend.title=element_blank(),
        axis.text.x = element_blank())+
  ggtitle("Females")+
  scale_colour_grey(start = 0, end = .45)


save(gwasACCres2, file = "results/GWAS/CO_count_QCed/Female-GWAS-lambda-corrected.RData")
save(gwasACC2, file = "results/GWAS/CO_count_QCed/Female-GWAS.RData")



# Female Only - Basic Model
sparrow.m <- sparrow %>% filter(SexF == "Male")

system.time(gwasACC_prefit2 <- preFitModel(fixed = yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2, 
                                           id.name = "id",
                                           genabel.data = data1,
                                           phenotype.data = sparrow.m, 
                                           corStruc = list(id = list("GRM")),
                                           GRM = sparrow.gkin.sym))

system.time(gwasACC2 <- rGLSadj(yapp_CO_count_QCed ~ 1 + Total_Coverage + Total_Coverage2,
                                genabel.data = data1,
                                phenotype.data = sparrow.m, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)
head(gwasACCres2)
plot_rGLSadj_results(gwasACCres2) +
  theme(axis.text = element_text(size = 24), 
        axis.title = element_text(size = 24,face = "bold"), 
        plot.title = element_text(size = 28),
        #plot.margin = margin(20, 50, 20, 50),
        legend.text=element_text(size=24),
        legend.title=element_blank(),
        axis.text.x = element_blank())+
  ggtitle("Males")+
  scale_colour_grey(start = 0, end = .45)


save(gwasACCres2, file = "results/GWAS/CO_count_QCed/Male-GWAS-lambda-corrected.RData")
save(gwasACC2, file = "results/GWAS/CO_count_QCed/Male-GWAS.RData")