# Calculate probability of having progressed disease under the state transition clock forward model

Calculates membership probability of the progressed disease state at a
given time (vectorized). This probability is from the state transition
clock forward model, according to the given statistical distributions
and parameters.

## Usage

``` r
prob_pd_stm_cf(time, dpam, starting = c(1, 0, 0))
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
prob_pd_stm_cf(0:100, params)
#>   [1] 0.000000000 0.057565843 0.102501319 0.140011833 0.171674410 0.198397032
#>   [7] 0.220826251 0.239465829 0.254684897 0.266784588 0.276064470 0.282811290
#>  [13] 0.287292879 0.289760924 0.290444507 0.289554760 0.287284771 0.283812230
#>  [19] 0.279299653 0.273905726 0.267784327 0.261061641 0.253856894 0.246277388
#>  [25] 0.238415524 0.230353913 0.222165759 0.213914700 0.205655880 0.197436766
#>  [31] 0.189297938 0.181273752 0.173393123 0.165680168 0.158154467 0.150831878
#>  [37] 0.143724753 0.136842555 0.130191950 0.123777527 0.117601707 0.111665238
#>  [43] 0.105967356 0.100506051 0.095278476 0.090280111 0.085506406 0.080951979
#>  [49] 0.076610924 0.072476967 0.068543515 0.064818251 0.061264396 0.057890339
#>  [55] 0.054688924 0.051652873 0.048764410 0.046038487 0.043456604 0.041012315
#>  [61] 0.038699151 0.036510856 0.034441387 0.032484908 0.030635803 0.028888675
#>  [67] 0.027238338 0.025679816 0.024208354 0.022819395 0.021508584 0.020271953
#>  [73] 0.019105157 0.018004616 0.016966732 0.015988084 0.015065419 0.014195646
#>  [79] 0.013376758 0.012604047 0.011875869 0.011189720 0.010543227 0.009934145
#>  [85] 0.009360347 0.008819822 0.008310665 0.007831078 0.007379362 0.006953909
#>  [91] 0.006553202 0.006175811 0.005820382 0.005485971 0.005170675 0.004873730
#>  [97] 0.004594067 0.004330694 0.004082633 0.003848996 0.003628937
# }
```
