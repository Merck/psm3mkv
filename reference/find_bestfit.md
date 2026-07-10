# Find the "best" survival regression from a list of model fits

When there are multiple survival regressions fitted to the same endpoint
and dataset, it is necessary to identify the preferred model. This
function reviews the fitted regressions and selects that with the
minimum Akaike or Bayesian Information Criterion (AIC, BIC), depending
on user choice. Model fits must be all parametric or all splines.

## Usage

``` r
find_bestfit(reglist, crit)
```

## Arguments

- reglist:

  List of fitted survival regressions to an endpoint and dataset.

- crit:

  Criterion to be used in selection of best fit, either "aic" (Akaike
  Information Criterion) or "bic" (Bayesian Information Criterion).

## Value

List of the single survival regression with the best fit.

## Examples

``` r
bosonc <- create_dummydata("flexbosms")
# Parametric modeling
fits_par <- fit_ends_mods_par(bosonc)
find_bestfit(fits_par$ttp, "aic")
#> $fit
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  -0.0376  -0.0652  -0.0101   0.0141
#> rate    0.0603   0.0455   0.0800   0.0087
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.3123, df = 2
#> AIC = 846.6246
#> 
#> 
#> $results
#> # A tibble: 6 × 13
#>      id valid conv  posdef  npts dists     pars loglik   aic   bic    ic rankaic
#>   <int> <lgl> <lgl> <lgl>  <dbl> <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1     1 TRUE  TRUE  TRUE     204 exp          1  -426.  853.  856.  853.       6
#> 2     2 TRUE  TRUE  TRUE     204 weibull…     2  -423.  850.  857.  850.       3
#> 3     3 TRUE  TRUE  TRUE     204 llogis       2  -422.  848.  855.  848.       2
#> 4     4 TRUE  TRUE  TRUE     204 lnorm        2  -424.  851.  858.  851.       5
#> 5     5 TRUE  TRUE  TRUE     204 gamma        2  -424.  851.  858.  851.       4
#> 6     6 TRUE  TRUE  TRUE     204 gompertz     2  -421.  847.  853.  847.       1
#> # ℹ 1 more variable: rankbic <dbl>
#> 
# \donttest{
# Splines modeling
fits_spl <- fit_ends_mods_spl(bosonc)
find_bestfit(fits_spl$ttp, "bic")
#> $fit
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -1.60314  -1.83224  -1.37404   0.11689
#> gamma1   0.38216   0.20690   0.55743   0.08942
#> gamma2  -0.01936  -0.03709  -0.00163   0.00905
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.5655, df = 3
#> AIC = 849.1309
#> 
#> 
#> $results
#> # A tibble: 9 × 14
#>      id valid conv  posdef  npts scales nknots  pars loglik   aic   bic    ic
#>   <int> <lgl> <lgl> <lgl>  <dbl> <chr>   <int> <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1     1 TRUE  TRUE  TRUE     204 hazard      1     3  -423.  852.  862.  862.
#> 2     2 TRUE  TRUE  TRUE     204 odds        1     3  -422.  849.  859.  859.
#> 3     3 TRUE  TRUE  TRUE     204 normal      1     3  -422.  849.  859.  859.
#> 4     4 TRUE  TRUE  TRUE     204 hazard      2     4  -423.  854.  867.  867.
#> 5     5 TRUE  TRUE  TRUE     204 odds        2     4  -421.  851.  864.  864.
#> 6     6 TRUE  TRUE  TRUE     204 normal      2     4  -421.  851.  864.  864.
#> 7     7 TRUE  TRUE  TRUE     204 hazard      3     5  -420.  850.  867.  867.
#> 8     8 TRUE  TRUE  TRUE     204 odds        3     5  -420.  850.  867.  867.
#> 9     9 TRUE  TRUE  TRUE     204 normal      3     5  -420.  850.  867.  867.
#> # ℹ 2 more variables: rankaic <dbl>, rankbic <dbl>
#> 
# }
```
