# Calculate probability of being alive under the state transition clock forward model

Calculates membership probability of being alive at a given time
(vectorized). This probability is from the state transition clock
forward model, according to the given statistical distributions and
parameters.

## Usage

``` r
prob_os_stm_cf(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized) from baseline.

- dpam:

  List of survival regressions for model endpoints. This must include
  pre-progression death (PPD), time to progression (TTP) and post
  progression survival calculated under the clock forward model
  (PPS-CF).

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
prob_os_stm_cf(0:100, params)
#>   [1] 1.00000000 0.99438187 0.98388723 0.97045962 0.95480239 0.93736449
#>   [7] 0.91847520 0.89839291 0.87733215 0.85547549 0.83298470 0.81000158
#>  [13] 0.78664786 0.76303528 0.73926205 0.71541805 0.69158476 0.66783743
#>  [19] 0.64424475 0.62087956 0.59781713 0.57511174 0.55281670 0.53097889
#>  [25] 0.50963525 0.48881683 0.46854870 0.44885041 0.42973671 0.41121802
#>  [31] 0.39330081 0.37598800 0.35927948 0.34317249 0.32766172 0.31273990
#>  [37] 0.29839789 0.28462513 0.27140962 0.25873850 0.24659790 0.23497329
#>  [43] 0.22384963 0.21321146 0.20304335 0.19332892 0.18405247 0.17519810
#>  [49] 0.16675002 0.15869261 0.15101046 0.14370296 0.13672548 0.13007897
#>  [55] 0.12374922 0.11772229 0.11197405 0.10651346 0.10131630 0.09637071
#>  [61] 0.09166506 0.08718819 0.08292942 0.07887846 0.07502548 0.07136107
#>  [67] 0.06787620 0.06456226 0.06141100 0.05841457 0.05556543 0.05285661
#>  [73] 0.05028089 0.04783192 0.04550350 0.04328969 0.04118485 0.03918361
#>  [79] 0.03728177 0.03547254 0.03375228 0.03211658 0.03056125 0.02908229
#>  [85] 0.02767592 0.02633852 0.02506668 0.02385712 0.02270676 0.02161266
#>  [91] 0.02057201 0.01958217 0.01864061 0.01774527 0.01689317 0.01608253
#>  [97] 0.01531128 0.01457750 0.01387931 0.01321495 0.01258276
# }
```
