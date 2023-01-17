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
# source("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/ASReml4.EstEffects.R")
source("C:/Users/johnb/Dropbox/McAuley PhD - Data/Scripts/Model/makeGRM.R")


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

# ___________________________________________________________________________
#  1. Setup eigenvalues
# ___________________________________________________________________________

# ___________________________________________________________________________
#  2. Perform Calculations
# ___________________________________________________________________________

my.n = #sample size of ind.
  my.m = #sample size of markers
  
  # Calculations effective sample sizze from eigenvalues of the kinship matrix
  # Assuming lambda = 1
  my.n0 = n0(lambda=1,eigvalues=d)

#effective correlation coefficient btw ind. in sample
rho(n=210,lambda=1,n0= my.n0)


power(n=210,h2=0.05,lambda=1,rho=,m=,alpha=0.05)


heritability1(n0=338.525,m=1619,lambda=1,power=0.85,alpha=0.05)


heritability2(n=210,rho=0.7651,m=1619,lambda=1,power=0.85,alpha=0.05)


sample(h2=0.05,rho=0.7651,m=1619,lambda=1,power=0.85,alpha=0.05)
