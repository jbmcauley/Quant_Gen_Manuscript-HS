## GWAS Power Analysis
## By John McAuley
## Jan 17th 2023


# ___________________________________________________________________________
#  Loading Libraries & Functions
# ___________________________________________________________________________

library(asreml)
library(tidyverse)
library(kinship2)
library(ggplot2)
library(reshape2)
library(cowplot)
library(glue)
library(dplyr)


n0<-function(lambda=1,eigvalues=d){
  return( (lambda+1)*sum(1/(eigvalues*lambda+1)) )
}

rho<-function(n=500,lambda=1,n0=500){
  y<-n0
  fn<-function(x,y){
    f<-(lambda+1)*( (n-1)/((1-x)*lambda+1) + 1/( (1+n*x-x)*lambda + 1)) - y
    return(f)
  }
  myrho<-uniroot(f=fn,y=y,lower=0,upper=1)
  return(myrho$root)
}

power<-function(n=500,h2=0.05,lambda=1,rho=0.5,m=1e5,alpha=0.05){
  alpha<-alpha/m
  x<-qchisq(1-alpha,1)
  n0<-(n-1)/((1-rho)*lambda+1) + 1/((1+n*rho-rho)*lambda+1)
  delta<-h2/(1-h2)*(lambda+1)*n0
  beta<-pchisq(x,1,delta)
  power<-1-beta
  return(power)
}

heritability1<-function(lambda=1,n0=2000,power=0.85,m=1e5,alpha=0.05){
  alpha<-alpha/m
  x<-qchisq(1-alpha,1)
  delta<-(qnorm(1-alpha/2)+qnorm(power))^2
  h2<-delta/(n0+delta)
  return(h2)
}


heritability2<-function(n=2000,lambda=1,rho=0.5,power=0.85,m=1e5,alpha=0.05){
  alpha<-alpha/m
  x<-qchisq(1-alpha,1)
  delta<-(qnorm(1-alpha/2)+qnorm(power))^2
  n0<-(lambda+1)*((n-1)/((1-rho)*lambda+1) + 1/((1+n*rho-rho)*lambda+1))
  h2<-delta/(n0+delta)
  return(h2)
}

sample<-function(h2=0.05,lambda=1,rho=0.5,power=0.85,m=1e5,alpha=0.05){
  alpha<-alpha/m
  y<-(qnorm(1-alpha/2)+qnorm(power))^2
  fn<-function(x,y){
    n0<-(x-1)/((1-rho)*lambda+1) + 1/((1+x*rho-rho)*lambda+1)
    f<-h2*(lambda+1)/(1-h2)*n0-y
    return(f)
  }
  myh2<-uniroot(f=fn,y=y,lower=0,upper=1e8)
  return(myh2$root)
}

# ___________________________________________________________________________
#  0. Load in Data & Setup variables
# ___________________________________________________________________________

# Make a kinship matrix

kk = read.csv(file = "data/GWAS/Kinship/kinship.csv", header = TRUE)
kk <- as.matrix(kk[,-1])
qq <- eigen(kk, symmetric = T)
d <- qq$values
u <- qq$vectors
write.csv(x = data.frame(d), "data/GWAS/Kinship/eigenvales.csv", row.names = F)
eigenvalues <- read.csv("data/GWAS/Kinship/eigenvales.csv", header = T)
d <- eigenvalues$d

# ___________________________________________________________________________
#  1. Perform Calculations
# ___________________________________________________________________________

my.n = 970    # sample size of ind.THIS VALUE DOES NOT ACCOUNT FOR LD btw neighboring markers
              # Need to estimate LD on dataset
my.m = 177909 # sample size of markers, autosomal markers from 180k snps
my.h2 = 0.05  # h2 of QTL
l = .1
# Calculations effective sample size from eigenvalues of the kinship matrix
# Assuming lambda = 1

n0(lambda=l,eigvalues=d)
my.n0 = n0(lambda=l,eigvalues=d)
# 1193.379

# effective correlation coefficient btw ind. in sample
rho(n=my.n,lambda=l,n0= my.n0)
my.rho <- rho(n=my.n,lambda=l,n0= my.n0)
# 0.376

#Power of sample n to detect QTL of h2 using m markers under rho
power(n = my.n*.5, h2= .05,lambda=l, rho = my.rho, m = my.m, alpha=0.05)

my.pwr <- power(n = my.n, h2= my.h2 ,lambda=l, rho = my.rho, m = my.m, alpha=0.05)

# Detectible h2 using effective sample size (n0)
heritability1(n0=my.n0,m=my.m,lambda=l,power=0.85, alpha=0.05)

# Detectible h2 using effective correlation (rho)
heritability2(n=my.n, rho=my.rho, m=my.m, lambda=l, power=0.85, alpha=0.05)

# minimum sample size required from this sample to detect a QTL with given h2 and power
sample(h2=my.h2, rho=my.rho, m=my.m, lambda=l, power=0.85, alpha=0.05)
