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
# RepeatABEL installation:
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





