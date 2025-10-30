library("data.table")
library("latticeExtra")
library(dplyr)
library(tidyverse)
library("zoo")

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

