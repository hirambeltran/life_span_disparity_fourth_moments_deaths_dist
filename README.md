# life_span_disparity_fourth_moments_deaths_dist
Computer code associated with the empirical estimation shown in the manuscript entitled "The Evolving Signature of Human Mortality: Skewness, Kurtosis, and Life Span Disparity"

There are three files that should be run sequentially, as indicated in the numbering of each file.
- Item 1 1_Edagger_discrete_and_trapezoids.r : computes edagger using trapezoids and a discrete aproximation, in addition to computing the first four moments of the deaths distribution for each country-year-sex in the dataset
- Item 2 2_Exact_derivatives_using_7th_degree_polynomials.r : fits a seventh degree polynomial to the life expectancy function and then estimates exact derivatives using coefficient estimates from the polynomial. This is done for each country-year-sex in the dataset
- Item 3 3_Decomposing_edagger.r : combines results from the previous two programs to empirically estimate the decomposition of edagger as a function of the four moments of the deaths distribution (see equation (2.5) in the paper). This program also empirically assesses the error term. Similarly, all calculations are done for each country-year-sex in the dataset.
