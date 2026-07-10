# Compare likelihoods of PSMs

Compare the total log-likelihood values for the patient-level dataset
after fitting PSM-simple and PSM-complex models to each combination of
endpoint distributions

## Usage

``` r
compare_psm_likes(ptdata, fitslist, cuttime = 0)
```

## Arguments

- ptdata:

  Dataset of patient level data. Must be a tibble with columns named:

  - ptid: patient identifier

  - pfs.durn: duration of PFS from baseline

  - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0
    for censoring)

  - os.durn: duration of OS from baseline

  - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)

  - ttp.durn: duration of TTP from baseline (usually should be equal to
    pfs.durn)

  - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for
    censoring).

- fitslist:

  List of distribution fits to relevant endpoints, after calling
  [`fit_ends_mods_par()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_par.md)
  or
  [`fit_ends_mods_spl()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_spl.md)

- cuttime:

  Time cutoff - this is nonzero for two-piece models.

## Value

List containing

- `results`: Dataset of calculation results for each model

- `bests`: Tibble indicating which is the best fitting model
  individually or jointly, to each endpoint, according to AIC or BIC

## Examples

``` r
# Fit parametric distributions to a dataset
bosonc <- create_dummydata("flexbosms")
parfits <- fit_ends_mods_par(bosonc)
# \donttest{
splfits <- fit_ends_mods_spl(bosonc)
# Present comparison of likelihood calculations
compare_psm_likes(bosonc, parfits)
#> Calculating PSM simple
#> Calculating PSM complex
#> Wrapping up
#> $results
#> # A tibble: 252 × 13
#>       id ttp_meth pfs_dist  os_dist      ll  npar  npts   aic   bic rank_aic
#>    <int> <chr>    <chr>     <chr>     <dbl> <dbl> <int> <dbl> <dbl>    <dbl>
#>  1     1 simple   exp       exp       -859.     3   204 1725. 1735.      139
#>  2     2 simple   exp       weibullPH -834.     4   203 1676. 1689.       38
#>  3     3 simple   exp       llogis    -839.     4   203 1687. 1700.       93
#>  4     4 simple   exp       lnorm     -839.     4   203 1686. 1699.       91
#>  5     5 simple   exp       gamma     -835.     4   203 1677. 1690.       42
#>  6     6 simple   exp       gompertz  -830.     4   203 1668. 1682.       24
#>  7     7 simple   weibullPH exp       -851.     4   204 1711. 1724.      126
#>  8     8 simple   weibullPH weibullPH -835.     5   203 1680. 1696.       53
#>  9     9 simple   weibullPH llogis    -841.     5   203 1692. 1709.      108
#> 10    10 simple   weibullPH lnorm     -839.     5   203 1687. 1704.       95
#> # ℹ 242 more rows
#> # ℹ 3 more variables: rank_bic <dbl>, best_aic <dbl>, best_bic <dbl>
#> 
#> $bests
#> # A tibble: 4 × 5
#>   ttp_meth pfs_dist os_dist   meth  ic   
#>   <chr>    <chr>    <chr>     <chr> <chr>
#> 1 gompertz exp      weibullPH ind   aic  
#> 2 gompertz exp      weibullPH ind   bic  
#> 3 gompertz lnorm    gamma     joint aic  
#> 4 gompertz lnorm    gamma     joint bic  
#> 
compare_psm_likes(bosonc, splfits)
#> Calculating PSM simple
#> Calculating PSM complex
#> Wrapping up
#> $results
#> # A tibble: 810 × 16
#>       id ttp_meth ttp_knots pfs_scales pfs_knots os_scales os_knots    ll  npar
#>    <int> <chr>        <dbl> <chr>          <int> <chr>        <int> <dbl> <dbl>
#>  1     1 simple           0 hazard             3 hazard           3 -835.     7
#>  2     2 simple           0 hazard             3 odds             3 -840.     7
#>  3     3 simple           0 hazard             3 normal           3 -838.     7
#>  4     4 simple           0 hazard             3 hazard           4 -833.     8
#>  5     5 simple           0 hazard             3 odds             4 -833.     8
#>  6     6 simple           0 hazard             3 normal           4 -832.     8
#>  7     7 simple           0 hazard             3 hazard           5 -833.     9
#>  8     8 simple           0 hazard             3 odds             5 -832.     9
#>  9     9 simple           0 hazard             3 normal           5 -832.     9
#> 10    10 simple           0 odds               3 hazard           3 -825.     7
#> # ℹ 800 more rows
#> # ℹ 7 more variables: npts <int>, aic <dbl>, bic <dbl>, rank_aic <dbl>,
#> #   rank_bic <dbl>, best_aic <dbl>, best_bic <dbl>
#> 
#> $bests
#> # A tibble: 4 × 8
#>   ttp_meth ttp_knots pfs_scales pfs_knots os_scales os_knots meth  ic   
#>   <chr>        <dbl> <chr>          <int> <chr>        <int> <chr> <chr>
#> 1 normal           3 normal             3 hazard           3 ind   aic  
#> 2 normal           3 normal             3 hazard           3 ind   bic  
#> 3 odds             3 normal             4 hazard           4 joint aic  
#> 4 simple           0 odds               3 normal           4 joint bic  
#> 
# }
```
