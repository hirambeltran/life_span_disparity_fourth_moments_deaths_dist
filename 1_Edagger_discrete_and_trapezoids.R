library("data.table")
library("latticeExtra")
library(stats)
library(foreign)
library(xtable)
library(stargazer)
library(dplyr)

library("purrr")
library(modelr)
library(tidyverse)

library("jtools")
library(dplyr)
library(broom)
library("huxtable")
library("zoo")

#setwd("C:/Users/Hiram/Box Sync/LAMBdA Dataset/Validation Tests LAMBdA_HMD_CDLT")



#load("Data/data_WHO_Hakkert_approx.RData")
#load("Data/data_HMD_Hakkert_approx.RData")
#load("Data/data_HMD_Hl_e0.RData")
#with(data.HakkertH,table(ctry))
#with(data.Hl.lne0,table(ctry))



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# 1. Human Mortality Database

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
setwd("E:/Box Sync/Beltran-Sanchez/Human Mort Database")
#setwd("C:/Users/beltrans/Box/Beltran-Sanchez/Human Mort Database")
load("HMD1x1_period_LT.RData")
hmd1x1<-data.table(hmd1x1.per)
hmd1x1[,list(min.yr=min(Year),max.yr=max(Year)),by=c("ctry")]

exclude<-c("ISL","FRACNP","DEUTNP","NZL_MA","NZL_NM","GBR_NP","GBRCENW","GBR_SCO","GBR_NIR")
hmd1x1<-subset(hmd1x1,!(ctry %in% exclude))
rm(exclude)
#View(hmd1x1[ctry=="FRATNP" & sex=="Female" & Year=="1816",])

# There is no life table data for Belgium in 1914-1919
hmd1x1<- hmd1x1[!(ctry=="BEL" & Year %in% c(1914:1919)),]
hmd1x1<- data.table(hmd1x1[!(ctry=="BEL" & Year %in% c(1914:1919)),])






setwd("E:/Box Sync/Beltran-Sanchez/Oscar Fernandez/Skewness-Kurtosis/")
#setwd("C:/Users/beltrans/Box/Beltran-Sanchez/Oscar Fernandez/Skewness-Kurtosis/")
library(zoo) # needed to compute rolling sums

save(hmd1x1 %>% filter(Year<=2019), file="Programs/ForRepository/HMD_period_upto2019.RData")
load()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Estimating e-dagger for all years, countries, sex using trapezoids

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
tmp<-lapply(0:100,FUN=function(a) {
  hmd1x1[Age>=a,] %>%   group_by(Year,ctry,sex) %>% 
         mutate(
         age = as.numeric(Age)) |>
         group_by(Year,ctry,sex) %>%   
         summarise(e.dagger.t = 0.5*(dx[1] * ex[1] + 
         		                       2*sum(dx[2:(length(Year)-1)] * ex[2:(length(Year)-1)]) + 
         		                       dx[length(Year)] * ex[length(Year)] ) / sum(dx)) %>%
        mutate(Age=a)
  })
edag.t<-data.table(do.call(rbind,tmp))


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Estimating e-dagger for all years, countries, sex using discrete approximation

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
tmp<-lapply(0:100,FUN=function(a) {
  hmd1x1[Age>=a,] %>%   group_by(Year,ctry,sex) %>% 
         mutate(
         age = as.numeric(Age),
         dx.ex=ifelse(Age==110,0,dx*0.5*rollapply(ex,2,sum))) %>%  
         group_by(Year,ctry,sex) %>%   summarise(e.dagger.d=sum(dx.ex)/sum(dx)) %>%
        mutate(Age=a)
  })
edag.d<-do.call(rbind,tmp)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Estimating sigma, mu3, and mu4 for ages 0 to 100 for all years, countries, sex
# using the trapezoidal rule

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
hmd1x1[,e_a:=ex+Age]
tmp<-lapply(0:100,FUN=function(a) {
  hmd1x1[Age>=a,] %>%  group_by(Year,ctry,sex) %>% 
        mutate(ex = ex[1],
  	     #age = ifelse(Age=="110+", "110", Age),
         age = as.numeric(Age),
         diff_sq = ((age - e_a[1])^2*dx),
         diff_3 = ((age - e_a[1])^3*dx),
         diff_4 = ((age - e_a[1])^4*dx)) %>% 
       group_by(Year,ctry,sex,ex) %>%   
        summarise(sigma2=0.5*(diff_sq[1]+2*sum(diff_sq[2:(length(Year)-1)])+diff_sq[length(Year)])/sum(dx),
                  mu3=0.5*(diff_3[1]+2*sum(diff_3[2:(length(Year)-1)])+diff_3[length(Year)])/sum(dx),
                  mu4=0.5*(diff_4[1]+2*sum(diff_4[2:(length(Year)-1)])+diff_4[length(Year)])/sum(dx)) %>%
        mutate(Age=a)
  })
sig.mu3.mu4<-data.table(do.call(rbind,tmp))



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# putting data all together

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
HMD.res<-merge(edag.t,edag.d, by=c("Year","ctry","sex","Age"))
HMD.res<-merge(HMD.res, sig.mu3.mu4, by=c("Year","ctry","sex","Age"))
HMD.res = data.table(setorder(HMD.res,ctry,sex,Year,Age))


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# computing mu3.hat and mu4.hat --> these are the standardized values

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
HMD.res[,sigma:=sqrt(sigma2)]
HMD.res[,mu3.hat:=mu3/(sigma^3)]
HMD.res[,mu4.hat:=mu4/(sigma^4)]
HMD.res

save(HMD.res,file="Data/HMD_var_edagger_by_ages_trapezoids.RData")


