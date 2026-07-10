# Calculate probability of being progression free in either state transition model (clock forward or clock reset)

Calculates membership probability for the progression free state, at a
particular time (vectorized), given either state transition model (clock
forward or clock reset) with given statistical distributions and
parameters.

## Usage

``` r
prob_pf_stm(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized)

- dpam:

  List of survival regressions for model endpoints. This must include
  pre-progression death (PPD) and time to progression (TTP).

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
prob_pf_stm(0:100, params)
#>   [1] 1.000000000 0.936816026 0.881385907 0.830447785 0.783127984 0.738967456
#>   [7] 0.697648951 0.658927085 0.622647253 0.588690906 0.556920232 0.527190295
#>  [13] 0.499354983 0.473274361 0.448817539 0.425863288 0.404299985 0.384025196
#>  [19] 0.364945100 0.346973830 0.330032807 0.314050098 0.298959805 0.284701505
#>  [25] 0.271219725 0.258462919 0.246382941 0.234935706 0.224080833 0.213781255
#>  [31] 0.204002873 0.194714251 0.185886357 0.177492325 0.169507255 0.161908026
#>  [37] 0.154673141 0.147782578 0.141217667 0.134960971 0.128996189 0.123308057
#>  [43] 0.117882272 0.112705414 0.107764877 0.103048810 0.098546060 0.094246124
#>  [49] 0.090139098 0.086215640 0.082466945 0.078884713 0.075461087 0.072188631
#>  [55] 0.069060301 0.066069417 0.063209645 0.060474974 0.057859696 0.055358390
#>  [61] 0.052965904 0.050677338 0.048488033 0.046393554 0.044389680 0.042472393
#>  [67] 0.040637862 0.038882440 0.037202649 0.035595173 0.034056849 0.032584661
#>  [73] 0.031175728 0.029827304 0.028536765 0.027301606 0.026119435 0.024987966
#>  [79] 0.023905015 0.022868496 0.021876415 0.020926865 0.020018024 0.019148149
#>  [85] 0.018315572 0.017518701 0.016756010 0.016026041 0.015327399 0.014658747
#>  [91] 0.014018808 0.013406359 0.012820230 0.012259301 0.011722498 0.011208797
#>  [97] 0.010717213 0.010246806 0.009796674 0.009365957 0.008953827
# }
```
