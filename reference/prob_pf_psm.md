# Calculate probability of being progression free in partitioned survival model

Calculates membership probability for the progression free state, at a
particular time (vectorized), given a partitioned survival model with
given statistical distributions and parameters.

## Usage

``` r
prob_pf_psm(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized)

- dpam:

  List of survival regressions for model endpoints. This must include
  progression-free survival (PFS).

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
prob_pf_psm(0:100, params)
#>   [1] 1.00000000 0.93912025 0.88389518 0.83206183 0.78346757 0.73798093
#>   [7] 0.69544066 0.65566985 0.61849121 0.58377788 0.55141951 0.52128082
#>  [13] 0.49321595 0.46707777 0.44272330 0.42001660 0.39883020 0.37904566
#>  [19] 0.36055364 0.34325356 0.32705314 0.31186789 0.29762051 0.28424030
#>  [25] 0.27166258 0.25982824 0.24868313 0.23817770 0.22826654 0.21890800
#>  [31] 0.21006386 0.20169901 0.19378114 0.18628054 0.17916982 0.17242372
#>  [37] 0.16601894 0.15993392 0.15414874 0.14864495 0.14340545 0.13841439
#>  [43] 0.13365703 0.12911971 0.12478968 0.12065511 0.11670495 0.11292892
#>  [49] 0.10931741 0.10586144 0.10255257 0.09938282 0.09634476 0.09343144
#>  [55] 0.09063634 0.08795338 0.08537683 0.08290135 0.08052190 0.07823376
#>  [61] 0.07603249 0.07391392 0.07187411 0.06990935 0.06801615 0.06619121
#>  [67] 0.06443143 0.06273386 0.06109573 0.05951441 0.05798741 0.05651238
#>  [73] 0.05508709 0.05370942 0.05237738 0.05108905 0.04984263 0.04863641
#>  [79] 0.04746875 0.04633811 0.04524300 0.04418203 0.04315386 0.04215722
#>  [85] 0.04119089 0.04025371 0.03934458 0.03846245 0.03760632 0.03677522
#>  [91] 0.03596824 0.03518449 0.03442315 0.03368341 0.03296451 0.03226570
#>  [97] 0.03158630 0.03092561 0.03028301 0.02965787 0.02904960
# }
```
