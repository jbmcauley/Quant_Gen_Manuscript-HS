## GWAS analysis for YAPP QC'd Dataset
## By John McAuley
## Dec 12th 2022
## R version 3.4.4

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
# Download package from CRAN archive
url <- "https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz"
pkgFile <- "GenABEL.data_1.0.0.tar.gz"
download.file(url = url, destfile = pkgFile)

# Install dependencies list in the DESCRIPTION file
# NA

# install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package
unlink(pkgFile)

# GenABEL
# Download package from CRAN archive
url <- "https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz"
pkgFile <- "GenABEL_1.8-0.tar.gz"
download.file(url = url, destfile = pkgFile)

#Look in the DESCRIPTION file to view dependencies needed.

# Install dependencies list in the DESCRIPTION file
install.packages(c("methods", "MASS", "utils"))
install.packages("qvalue", "DatABEL" , "hglm",
"MetABEL", "PredictABEL", "VariABEL", "bigRR")
# install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package
unlink(pkgFile)


#_______________________________________________________
# 1.  Setup 
#_______________________________________________________
# Requirements: Phenotype Data, Genetic data, and GRM
# load rGLSadj.R as it contains functions we will use
source("rGLSadj.R")

# Load in a dataframe with Phenotypic data
load("...")

# Load in genotype data
load("...")

# gkin setup
sparrow.gkin.sym <- sparrow.gkin
sparrow.gkin.sym[upper.tri(sparrow.gkin.sym)] = t(sparrow.gkin.sym)[upper.tri(sparrow.gkin.sym)]
sparrow.gkin.sym <- sparrow.gkin.sym * 2

#Inspect genotype object
qc0snp <- summary.snp.data(sparrowgen@gtdata)
hist(qc0snp$CallRate, breaks = 50)
qc0id <- perid.summary(sparrowgen)
hist(qc0id$CallPP, breaks = 50) 

#Choose Quality Control thresholds for check.marker() function based on histograms
qc1 <- check.marker(sparrowgen, callrate = 0.9, perid.call = 0.8, p.level = 0)
data1 <- sparrowgen[qc1$idok, qc1$snpok]
rm(sparrowgen)
rm(qc1)


#_______________________________________________________
# 2.  RepeatABEL 
#_______________________________________________________
system.time(gwasACC_prefit2 <- preFitModel(fixed = ACC ~ sex, 
                                           id.name = "id",
                                           genabel.data = data1,
                                           phenotype.data = sparrow, 
                                           corStruc = list(id = list("GRM")),
                                           GRM = sparrow.gkin.sym))
system.time(gwasACC2 <- rGLSadj(ACC ~ sex,
                                genabel.data = data1,
                                phenotype.data = sparrow, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)
head(gwasACCres2)
plot_rGLSadj_results(test) +
  theme(axis.text = element_text(size = 24), 
        axis.title = element_text(size = 24,face = "bold"), 
        plot.title = element_text(size = 28),
        #plot.margin = margin(20, 50, 20, 50),
        legend.text=element_text(size=24),
        legend.title=element_blank(),
        axis.text.x = element_blank())+
  ggtitle("Both Sexes")+
  scale_colour_grey(start = 0, end = .45)


save(gwasACCres2, file = "All-sparrow-GWAS-lambda-corrected.RData")
save(gwasACCres2, file = "All-sparrow-GWAS.RData")
