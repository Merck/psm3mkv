# Calculate probability of post progression survival under the state transition clock reset model

Calculates probability of post progression survival at a given time from
progression (vectorized). This probability is from the state transition
clock reset model, according to the given statistical distributions and
parameters.

## Usage

``` r
prob_pps_cr(time, dpam)
```

## Arguments

- time:

  Time (numeric and vectorized) from baseline - not time from
  progression.

- dpam:

  List of survival regressions for model endpoints. This must include
  post progression survival calculated under the clock reset state
  transition model.

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
prob_pps_cr(0:100, params)
#>   [1] 1.000000000 0.987105716 0.934461106 0.868401686 0.806475394 0.753434376
#>   [7] 0.708946460 0.669965664 0.634136012 0.600110241 0.567138271 0.534824586
#>  [13] 0.502988745 0.471688904 0.441242337 0.411907557 0.383860425 0.357212477
#>  [19] 0.332025219 0.308321522 0.286094792 0.265316380 0.245941526 0.227914131
#>  [25] 0.211170533 0.195642481 0.181259439 0.167950356 0.155644991 0.144274888
#>  [31] 0.133774071 0.124079519 0.115131461 0.106873544 0.099252896 0.092220108
#>  [37] 0.085729164 0.079737331 0.074205016 0.069095609 0.064375317 0.060012987
#>  [43] 0.055979934 0.052249773 0.048798245 0.045602808 0.042642502 0.039898200
#>  [49] 0.037352466 0.034989399 0.032794482 0.030754457 0.028857209 0.027091656
#>  [55] 0.025447658 0.023915928 0.022487958 0.021155944 0.019912726 0.018751730
#>  [61] 0.017666915 0.016652728 0.015704061 0.014816209 0.013984840 0.013205960
#>  [67] 0.012475887 0.011791221 0.011148823 0.010545793 0.009979448 0.009447308
#>  [73] 0.008947078 0.008476629 0.008033991 0.007617337 0.007224971 0.006855320
#>  [79] 0.006506923 0.006178421 0.005868551 0.005576140 0.005300095 0.005039397
#>  [85] 0.004793098 0.004560316 0.004340227 0.004132060 0.003935100 0.003748675
#>  [91] 0.003572161 0.003404972 0.003246561 0.003096416 0.002954059 0.002819042
#>  [97] 0.002690943 0.002569370 0.002453954 0.002344348 0.002240229
# }
```
