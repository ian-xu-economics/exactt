
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
#> Package 'exactt' | Version 1.2.7
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

print(exactt.1, digits = 5)
#> 
#> Call:
#> exactt(model = len ~ dose + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1)
#> 
#> 
#> Summary:
#>         Estimate  Pr(>|t|)     5% W   95% W       5%     95%
#> dose      9.7636   0.07500    2.414  16.517    2.414  16.517
#> suppVC   -3.7000   0.26667  -11.880  10.300  -11.880  10.300
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

print(exactt.1, digits = 5)
#> 
#> Call:
#> exactt(model = len ~ dose + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1)
#> 
#> 
#> Summary:
#>         Estimate  Pr(>|t|)     5% W   95% W       5%     95%
#> dose      9.7636   0.07500    2.414  16.517    2.414  16.517
#> suppVC   -3.7000   0.26667  -11.880  10.300  -11.880  10.300
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
#>                   Estimate  Pr(>|t|)  5% W  95% W    5%   95%
#> as.factor(dose)1      9.13      0.05  5.44   70.3  5.44  70.3
#> as.factor(dose)2     15.49      1.00  -Inf    Inf  -Inf   Inf
#> suppVC               -3.70      1.00  -Inf    Inf  -Inf   Inf
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

print(exactt.4, digits = 5)

#> ✔ Optimizing ordering for `as.factor(dose)1`.
#> GA | iter = 1 | Mean = 2.992124 | Best = 4.942080
#> GA | iter = 2 | Mean = 3.134706 | Best = 4.942080
#> GA | iter = 3 | Mean = 3.144138 | Best = 5.381865
#> GA | iter = 4 | Mean = 3.069222 | Best = 5.381865
#> GA | iter = 5 | Mean = 2.897159 | Best = 5.381865
#> ✔ Optimizing ordering for `as.factor(dose)2`.
#> GA | iter = 1 | Mean = 2.799936 | Best = 5.007930
#> GA | iter = 2 | Mean = 3.107845 | Best = 5.244658
#> GA | iter = 3 | Mean = 3.469066 | Best = 5.244658
#> GA | iter = 4 | Mean = 3.356677 | Best = 5.244658
#> GA | iter = 5 | Mean = 3.129926 | Best = 5.244658
#> ✔ Optimizing ordering for `suppVC`.
#> GA | iter = 1 | Mean = 4.164565 | Best = 6.068181
#> GA | iter = 2 | Mean = 4.608127 | Best = 7.120173
#> GA | iter = 3 | Mean = 4.861197 | Best = 7.120173
#> GA | iter = 4 | Mean = 4.465345 | Best = 7.791005
#> GA | iter = 5 | Mean = 4.438026 | Best = 7.791005
#> 
#> Call:
#> exactt(model = len ~ as.factor(dose) + supp, data = datasets::ToothGrowth, 
#>     alpha = 0.1, optimize = TRUE, seed = 2024, parallel = TRUE, 
#>     maxiter = 5)
#> 
#> 
#> Summary:
#>                   Estimate   Pr(>|t|)    5% W   95% W      5%     95%
#> as.factor(dose)1     9.130  0.0083333   5.120  14.600   5.120  14.600
#> as.factor(dose)2    15.495  0.0083333  12.675  19.750  12.670  19.750
#> suppVC              -3.700  0.0583330  -7.044  -1.085  -7.044  -1.085
```

Note that by optimizing the data ordering, `exactt()` is now able to
construct informative 90% confidence intervals for each category of
`dose` and `supp` when they equal “2” and “VC” respectively.
Furthermore, the detailed results of the optimization process, including
the genetic algorithm’s configurations and outcomes for each variable,
are stored in the `exactt.4$gaResults`. For instance, to review a
summary of the genetic algorithm’s performance for the `suppVC`
variable, use:

``` r
exactt.4$gaResults$suppVC@summary

#>           max     mean       q3   median       q1      min
#> [1,] 6.068181 4.164565 5.062314 4.052930 3.469694 2.298845
#> [2,] 7.120173 4.608127 5.445378 4.569416 3.743193 1.569542
#> [3,] 7.120173 4.861197 5.918520 5.205500 3.757399 1.944108
#> [4,] 7.791005 4.465345 5.342829 4.491921 3.364883 1.883697
#> [5,] 7.791005 4.438026 5.267938 4.423997 3.567307 1.555700
```

### Note on Optimization Effects

While optimization generally improves statistical power, it is essential
to remember that it increases the average power and may not universally
reduce the confidence interval’s width in every instance.

## Example Usage: IV Case

The `exactt()` function is capable of handling models with instrumental
variables (IV). In Example 15.5 of Wooldridge (2020), Wooldridge
reanalyzes Mroz (1987). This example explores the impact of education
(`educ`) on `log(wage)`, using parental education levels—mother’s
education (`motheduc`) and father’s education (`fatheduc`)—as
instruments. The model controls for experience (`exper`) and its square
(`expersq`), with education being the primary variable of interest,
hence we set variables = 1. Optionally, as before, we can optimize the
data ordering to enhance statistical power.

``` r
exactt.iv <- exactt(lwage ~ educ + exper + expersq | exper + expersq + motheduc + fatheduc,
                    data = wooldridge::mroz,
                    variables = 1,
                    optimize = TRUE,
                    parallel = TRUE,
                    maxiter = 10,
                    monitor = TRUE,
                    seed = 2024,
                    GX1 = FALSE,
                    studentize = FALSE)

exactt.iv

#> ✔ Optimizing ordering for `educ`.
#> GA | iter = 1 | Mean = 2717502 | Best = 3126138
#> GA | iter = 2 | Mean = 2771190 | Best = 3214122
#> GA | iter = 3 | Mean = 2836298 | Best = 3214122
#> GA | iter = 4 | Mean = 2793909 | Best = 3251129
#> GA | iter = 5 | Mean = 2822831 | Best = 3251129
#> GA | iter = 6 | Mean = 2789954 | Best = 3320237
#> GA | iter = 7 | Mean = 2798405 | Best = 3320237
#> GA | iter = 8 | Mean = 2819648 | Best = 3320237
#> GA | iter = 9 | Mean = 2797339 | Best = 3320237
#> GA | iter = 10 | Mean = 2794273 | Best = 3320237
#> 
#> Call:
#> exactt(model = lwage ~ educ + exper + expersq | exper + expersq + 
#>     motheduc + fatheduc, data = wooldridge::mroz, variables = 1, 
#>     studentize = FALSE, optimize = TRUE, GX1 = FALSE, seed = 2024, 
#>    parallel = TRUE, maxiter = 10, monitor = TRUE)
#> 
#> 
#> Summary:
#>       Estimate  Pr(>|t|)   2.5% W  97.5% W     2.5%  97.5%
#> educ    0.0614   0.06667  -0.0057    0.152  -0.0057  0.152
```
