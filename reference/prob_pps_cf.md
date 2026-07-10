# Calculate probability of post progression survival under the state transition clock forward model

Calculates probability of post progression survival at a given time from
progression (vectorized). This probability is from the state transition
clock forward model, according to the given statistical distributions
and parameters.

## Usage

``` r
prob_pps_cf(ttptimes, ppstimes, dpam)
```

## Arguments

- ttptimes:

  Time (numeric and vectorized) from progression - not time from
  baseline.

- ppstimes:

  Time (numeric and vectorized) of progression

- dpam:

  List of survival regressions for model endpoints. This must include
  post progression survival calculated under the clock forward state
  transition model.

## Value

Vector of the mean probabilities of post-progression survival at each
PPS time, averaged over TTP times.

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
prob_pps_cf(0:100, 0:100, params)
#>   [1] 1.0000000000 0.9269277329 0.8595239740 0.7969376513 0.7387652556
#>   [6] 0.6846798786 0.6343928456 0.5876418128 0.5441853505 0.5037998452
#>  [11] 0.4662774338 0.4314244718 0.3990603131 0.3690162875 0.3411348163
#>  [16] 0.3152686285 0.2912800553 0.2690403877 0.2484292861 0.2293342129
#>  [21] 0.2116498123 0.1952774482 0.1801248358 0.1661056875 0.1531393715
#>  [26] 0.1411505858 0.1300690452 0.1198291833 0.1103698686 0.1016341336
#>  [31] 0.0935689188 0.0861248289 0.0792559018 0.0729193903 0.0670755557
#>  [36] 0.0616874726 0.0567208450 0.0521438326 0.0479268878 0.0440426013
#>  [41] 0.0404655582 0.0371722017 0.0341407061 0.0313508565 0.0287839376
#>  [46] 0.0264226280 0.0242509024 0.0222539395 0.0204180363 0.0187305277
#>  [51] 0.0171797114 0.0157547788 0.0144457490 0.0132434086 0.0121392552
#>  [56] 0.0111254450 0.0101947448 0.0093404857 0.0085565206 0.0078371848
#>  [61] 0.0071772589 0.0065719348 0.0060167848 0.0055077315 0.0050410215
#>  [66] 0.0046132003 0.0042210887 0.0038617618 0.0035325292 0.0032309164
#>  [71] 0.0029546483 0.0027016332 0.0024699484 0.0022578272 0.0020636461
#>  [76] 0.0018859137 0.0017232603 0.0015744280 0.0014382617 0.0013137012
#>  [81] 0.0011997733 0.0010955850 0.0010003168 0.0009132171 0.0008335964
#>  [86] 0.0007608223 0.0006943151 0.0006335432 0.0005780195 0.0005272975
#>  [91] 0.0004809679 0.0004386560 0.0004000183 0.0003647404 0.0003325341
#>  [96] 0.0003031359 0.0002763042 0.0002518181 0.0002294753 0.0002090907
#> [101] 0.0001904951
# }
```
