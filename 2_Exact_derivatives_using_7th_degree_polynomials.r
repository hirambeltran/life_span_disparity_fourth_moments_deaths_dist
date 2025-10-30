library("data.table")
library(dplyr)
library(Hmisc)
library(tidyverse)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# 1. Human Mortality Database

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
hmd1x1[ , e_a := ex + Age]


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Estimate a 7-th degree polynomial to e(x)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
poly7<-hmd1x1 %>% nest_by(Year,ctry,sex) %>% 
	       mutate(pol7 = list(lm( ex ~ 1 + poly(Age, 7, raw=TRUE), data=data)))




# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Evaluating each polynomial at value e_a

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
poly7.e_a=data.table(poly7 %>% 
		       reframe(Age = data$Age, ex.e_a = predict(pol7, newdata = data.frame(Age = data$e_a))) )

# poly7.e_a[ctry=="FRATNP" & sex=="Female" & Year=="1816" ,]
#anyDuplicated(poly7.e_a[, .(ctry, Year, sex, Age)])
#poly7.e_a<- unique(poly7.e_a)

# putting the data together
hmd1x1 <- left_join(hmd1x1, poly7.e_a, by = c("Year","ctry","sex","Age"), keep=FALSE)





# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Computing EXACT derivatives of the polynomial evaluated at value e_a

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# extracting the coeff from each polynomial estimate
t1<- data.table(poly7 %>% reframe(tidy(pol7)) )

# coefficients associated with "x" for the derivatives of the 7-th degree polynomial
# f(x)=   a + b1x + b2x^2 + b3x^3   + b4x^4    + b5x^5    + b6x^6    + b7x^7
# f'(1) = a + b1  + 2*b2x + 3*b3x^2 + 4*b4x^3  + 5*b5x^4  + 6*b6x^5  + 7*b7x^6
# f''(1) =           2*b2 + 6*b3x   + 12*b4x^2 + 20*b5x^3 + 30*b6x^4 + 42*b7x^5 --> coef: c(0,0,2,6,12,20,30,42)
t1[,d2.coef:= c(0,0,2,6,12,20,30,42), by=c("Year","ctry","sex")]

# f3(1) =                    6*b3   + 24*b4x   + 60*b5x^2 + 120*b6x^3 + 210*b7x^4 --> coef: c(0,0,0,6,24,60,120,210)
t1[,d3.coef:= c(0,0,0,6,24,60,120,210), by=c("Year","ctry","sex")]

# f4(1) =                           + 24*b4    + 120*b5x  + 360*b6x^2 + 840*b7x^3 --> coef: c(0,0,0,0,24,120,360,840)
t1[,d4.coef:= c(0,0,0,0,24,120,360,840), by=c("Year","ctry","sex")]
t1[,id:=paste0(Year,ctry,sex)]

# actual computation of the EXACT derivatives
res<-lapply(unique(t1$id),FUN=function(s){
  # select a ctry-year-sex, example, s="1816FRATNPFemale"
  coeffs<-t1[id==s,]
  hmd<-hmd1x1[id==s,]
  
  # e_a raised to the power of each x in the corresponding derivative order
  derivs<-lapply(hmd$e_a,FUN=function(x){  
    data.table(e_a=x,
               d2.e_a=sum(coeffs$estimate*coeffs$d2.coef*c(0,0,1,x,x^2,x^3,x^4,x^5)),
               d3.e_a=sum(coeffs$estimate*coeffs$d3.coef*c(0,0,0,1,x,x^2,x^3,x^4)),
               d4.e_a=sum(coeffs$estimate*coeffs$d4.coef*c(0,0,0,0,1,x,x^2,x^3)))
  })  
  out<- do.call(rbind,derivs) %>% mutate(Age=hmd$Age,id=s)
})  
poly7.deriv.e_a <- data.table(do.call(rbind,res))




# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Putting the results together: 
#    evaluated polynomial value of e_a and EXACT derivatives of the 7th degree polynomial also evaluaated at e_a

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
poly7.deriv.e_a <- left_join(hmd1x1, poly7.deriv.e_a[ , .(id, Age, d2.e_a, d3.e_a, d4.e_a)], by=c("id","Age")) 



