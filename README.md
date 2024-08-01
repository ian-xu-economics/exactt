
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exactt <img src="man/figures/package-sticker.png" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/exactt)](https://CRAN.R-project.org/package=exactt)
[![R CMD
Check](https://github.com/ian-xu-economics/exactt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ian-xu-economics/exactt/actions/workflows/R-CMD-check.yaml)
[![Website](https://github.com/ian-xu-economics/exactt/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ian-xu-economics/exactt/actions/workflows/pkgdown.yaml)
[![Test
coverage](https://github.com/ian-xu-economics/exactt/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/ian-xu-economics/exactt/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/ian-xu-economics/exactt/branch/main/graph/badge.svg)](https://codecov.io/gh/ian-xu-economics/exactt)
<!-- badges: end -->

## Introduction

The `exactt` package tests whether a slope coefficient is equal to some
null value using the novel method described in Pouliot (2023).
Importantly, inverting such a test produces a marginally valid
confidence interval.

## Installation

The `exactt` package is hosted on GitHub at
<https://github.com/ian-xu-economics/exactt/>. It can be installed using
the `remotes::install_github()` function:

``` r
# install.packages("remotes")
remotes::install_github("ian-xu-economics/exactt")
```

## Attribution

To cite the `exactt` package in publications, use the `citation()`
function, which provides both the text version and the BibTeX entry for
referencing:

``` r
citation("exactt")
```

## Using `exactt`

After installing `exactt`, we can attach the package to our session
using the base `library()` function:

``` r
library("exactt")
#> Package 'exactt' | Version 1.2.1
```

## Example Usage: Regular Case

To compute the $(1-\alpha)$ confidence interval, use the `exactt()`
function. Here’s an example looking at the effect of vitamin C on tooth
growth in guinea pigs using data from `datasets::ToothGrowth`. We’ll
investigate the relationship between `supp` (orange juice (OJ) or
ascorbic acid (VC)) and `dose` (dose in milligrams/day) on `len` (tooth
length).

``` r
summary(datasets::ToothGrowth)
#>       len        supp         dose      
#>  Min.   : 4.20   OJ:30   Min.   :0.500  
#>  1st Qu.:13.07   VC:30   1st Qu.:0.500  
#>  Median :19.25           Median :1.000  
#>  Mean   :18.81           Mean   :1.167  
#>  3rd Qu.:25.27           3rd Qu.:2.000  
#>  Max.   :33.90           Max.   :2.000
```

Suppose our model is
$len_i = \beta_0 + \beta_{dose} \times dose_i + \beta_{supp} \times supp_i + \varepsilon_i$.
We can create a 90% confidence interval by plugging in standard formula
notation into `exactt()`. The level of significance (alpha) equals 0.1
here, but if we choose not to specify any additional parameters, then by
default:

- The number of blocks used equals 5 (`nBlocks = 5`).
- The confidence interval is constructed for all variables
  (`variables = NULL`).
- The number of permutations is equal to
  (`nPerms = factorial(nBlocks)`).
- The level of significance equals 0.05 (`alpha = 0.05`).
- The test statistics are studentized (`studentize = TRUE`).
- The ordering of the data is not permuted (`permutation = NULL`).
- The ordering of the data is not optimized (`optimize = FALSE`).

``` r
exactt.1 <- exactt(model = len ~ dose + supp,
                   data = datasets::ToothGrowth,
                   alpha = 0.1)

exactt.1
#> 
#> Call:
#> exactt(model = len ~ dose + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1)
#> 
#> 
#> Summary:
#>         Estimate  Pr(>|t|)     5% W   95% W       5%    95%
#> dose       9.764    0.0750    2.414  16.520    2.414  16.52
#> suppVC    -3.700    0.2583  -10.980   7.608  -10.990   7.61
```

## Focusing on Specific Variables

To focus on specific coefficients, set the `variables` parameter. The
number entered corresponds to the index of the regressors in the model
(note that the intercept is never counted). For example, set
`variables = 1` for `dose`, and set `variables = 2` for `supp`.

``` r
exactt.2 <- exactt(model = len ~ dose + supp,
                   data = datasets::ToothGrowth,
                   alpha = 0.1,
                   variables = 1)

exactt.2
#> 
#> Call:
#> exactt(model = len ~ dose + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1, variables = 1)
#> 
#> 
#> Summary:
#>       Estimate  Pr(>|t|)   5% W  95% W     5%    95%
#> dose     9.764     0.075  2.414  16.52  2.414  16.52
```

This creates a 90% confidence interval for `dose` only. It is equivalent
to the case where `variables = NULL` (all variables are of interest)
because these confidence intervals are marginally valid.

## Model Flexibility

The `exactt()` function is designed to allow for easy modification of
your model. For instance, you can treat a variable as categorical,
include polynomial terms, or apply other transformations directly within
the model formula. This flexibility helps tailor the analysis to
specific research questions without needing pre-transformed data. To
illustrate, consider treating `dose` as a categorical variable to
explore its discrete impact on tooth length:

``` r
exactt.3 <- exactt(model = len ~ as.factor(dose) + supp,
                   data = datasets::ToothGrowth,
                   alpha = 0.1)

exactt.3
#> 
#> Call:
#> exactt(model = len ~ as.factor(dose) + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1)
#> 
#> 
#> Summary:
#>                   Estimate  Pr(>|t|)   5% W  95% W    5%    95%
#> as.factor(dose)1      9.13   0.01667  5.647  19.13  5.64  19.13
#> as.factor(dose)2     15.49   0.05833  9.471    Inf  9.47    Inf
#> suppVC               -3.70   0.23330   -Inf    Inf  -Inf    Inf
```

The 90% confidence intervals when `dose` equals “2” and `supp` equals
“VC” is not informative due to suboptimal data ordering, which can
diminish the statistical power of the test. This issue can be addressed
by optimizing the data ordering.

## Optimizing Data Ordering

The confidence intervals produced by the `exactt()` function can change
with the ordering of the data. Certain data orderings can enhance
statistical power, particularly when the sample size is small and the
number of blocks is large. The impact of optimization is even more
pronounced when dealing with categorical variables, where appropriate
ordering can substantially increase the test’s power.

The `exactt()` function utilizes a genetic algorithm (provided by the
`GA::ga()` function) to optimize data ordering. This approach
systematically explores various data arrangements to find the one that
maximizes statistical power on average.

## Enabling Optimization

To activate the optimization feature, set `optimize = TRUE`.
Additionally, `exactt()` allows for the specification of various
parameters of the `GA::ga()` function to tailor the optimization
process. For instance, you can limit the number of iterations with
`maxiter` or specify the seed with `seed` for reproducibility:

``` r
exactt.4 <- exactt(model = len ~ as.factor(dose) + supp,
                   data = datasets::ToothGrowth,
                   alpha = 0.1,
                   optimize = TRUE,
                   parallel = TRUE,
                   maxiter = 5,
                   seed = 2024)

exactt.4

#> GA | iter = 1 | Mean = 2.992124 | Best = 4.942080
#> GA | iter = 2 | Mean = 3.109950 | Best = 5.195252
#> GA | iter = 3 | Mean = 2.982926 | Best = 5.195252
#> GA | iter = 4 | Mean = 3.037489 | Best = 5.195252
#> GA | iter = 5 | Mean = 3.053269 | Best = 5.418928
#> GA | iter = 1 | Mean = 2.799936 | Best = 5.007930
#> GA | iter = 2 | Mean = 3.017456 | Best = 5.007930
#> GA | iter = 3 | Mean = 3.140062 | Best = 5.007930
#> GA | iter = 4 | Mean = 3.105836 | Best = 5.007930
#> GA | iter = 5 | Mean = 3.191111 | Best = 5.007930
#> GA | iter = 1 | Mean = 4.164565 | Best = 6.068181
#> GA | iter = 2 | Mean = 4.263754 | Best = 6.343909
#> GA | iter = 3 | Mean = 4.649025 | Best = 7.158755
#> GA | iter = 4 | Mean = 4.668351 | Best = 7.158755
#> GA | iter = 5 | Mean = 4.520408 | Best = 7.158755
#> 
#> Call:
#> exactt(model = len ~ as.factor(dose) + supp, data = datasets::ToothGrowth, 
#>     optimize = TRUE, parallel = TRUE, maxiter = 5, seed = 2024)
#> 
#> 
#> Summary:
#>                   Estimate  Pr(>|t|)    5% W   95% W      5%     95%
#> as.factor(dose)1      9.13  0.008333   8.030  11.420   8.030  11.420
#> as.factor(dose)2     15.49  0.016670  12.160  19.920  12.160  19.920
#> suppVC               -3.70  0.033330  -6.687  -1.246  -6.687  -1.246
```
