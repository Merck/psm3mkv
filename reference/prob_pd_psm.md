# Calculate membership probability of progressed disease state in a partitioned survival model

Calculates membership probability of having progressed disease at a
particular time (vectorized), given the partitioned survival model with
certain statistical distributions and parameters.

## Usage

``` r
prob_pd_psm(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized)

- dpam:

  List of survival regressions for model endpoints. This must include
  progression-free survival (PFS) and overall survival (OS).

- starting:

  Vector of membership probabilities (PF, PD, death) at time zero.

## Value

Numeric value

## Examples

``` r
# \donttest{
bosonc <- create_dummydata("flexbosms")
fits <- fit_ends_mods_spl(bosonc)
# Pick out best distribution according to min AIC
params <- list(
  ppd = find_bestfit(fits$ppd, "aic")$fit,
  ttp = find_bestfit(fits$ttp, "aic")$fit,
  pfs = find_bestfit(fits$pfs, "aic")$fit,
  os = find_bestfit(fits$os, "aic")$fit,
  pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
  pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
)
prob_pd_psm(0:100, params)
#>   [1]  0.000000000  0.054213562  0.098565251  0.137013060  0.170266758
#>   [6]  0.198796057  0.223008955  0.243280274  0.259955242  0.273308351
#>  [11]  0.283582761  0.291034994  0.295922768  0.298496989  0.298997351
#>  [16]  0.297650474  0.294672507  0.290266681  0.284620862  0.277908013
#>  [21]  0.270286856  0.261902614  0.252887788  0.243362937  0.233437439
#>  [26]  0.223210226  0.212770498  0.202198385  0.191565592  0.180935996
#>  [31]  0.170366221  0.159906170  0.149599537  0.139484288  0.129593106
#>  [36]  0.119953824  0.110589822  0.101520409  0.092761177  0.084324333
#>  [41]  0.076219022  0.068451612  0.061025979  0.053943759  0.047204597
#>  [46]  0.040806365  0.034745377  0.029016580  0.023613739  0.018529600
#>  [51]  0.013756110  0.009284544  0.005105561  0.001209324 -0.002414425
#>  [56] -0.005776271 -0.008886978 -0.011757406 -0.014398443 -0.016820940
#>  [61] -0.019035660 -0.021053227 -0.022884090 -0.024538479 -0.026026384
#>  [66] -0.027357522 -0.028541320 -0.029586897 -0.030503052 -0.031298251
#>  [71] -0.031980622 -0.032557951 -0.033037680 -0.033426907 -0.033732387
#>  [76] -0.033960539 -0.034117447 -0.034208871 -0.034240252 -0.034216719
#>  [81] -0.034143102 -0.034023939 -0.033863488 -0.033665734 -0.033434406
#>  [86] -0.033172981 -0.032884702 -0.032572583 -0.032239423 -0.031887818
#>  [91] -0.031520167 -0.031138688 -0.030745425 -0.030342259 -0.029930916
#>  [96] -0.029512981 -0.029089899 -0.028662992 -0.028233462 -0.027802399
#> [101] -0.027370792
# }
```
