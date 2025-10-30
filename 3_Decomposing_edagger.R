library("data.table")
library(dplyr)
library(tidyverse)



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Decomposing e-dagger into 4 terms plus a residual 
# edagger(a)= e(e_a) + (d2.e_a/2!)*sigma^2 + (d3.e_a/3!)*sigma^3*skewness + (d4.e_a/4!)*sigma^4*kurtosis + R_a(a)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# file created in R-code "1_Edagger_discrete_and_trapezoids_and_4moments.r"
HMD.res

# file created in R-code "2_Exact_derivatives_using_7th_degree_polynomials.r"
poly7.deriv.e_a        

all<- left_join(HMD.res, 
	              poly7.deriv.e_a[ , .(Year,ctry,sex,Age,e_a,ex.e_a,d2.e_a,d3.e_a,d4.e_a)], 
	              by=c("Year","ctry","sex","Age"))

all[ , `:=`  (term1=ex.e_a,
              term2=(d2.e_a/factorial(2)) * sigma^2,
              term3=(d3.e_a/factorial(3)) * (sigma^3) * mu3.hat,
              term4=(d4.e_a/factorial(4)) * (sigma^4) * mu4.hat)]

all[ , terms.sum := term1 + term2 + term3 + term4]
all[ , error := e.dagger - (term1 + term2 + term3 + term4), by=c("Year","ctry","sex")]
all[ , rel.error := error/e.dagger, by = c("Year","ctry","sex")]


 