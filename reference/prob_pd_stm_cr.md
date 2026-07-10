# Calculate probability of having progressed disease under the state transition clock reset model

Calculates membership probability of the progressed disease state at a
given time (vectorized). This probability is from the state transition
clock reset model, according to the given statistical distributions and
parameters.

## Usage

``` r
prob_pd_stm_cr(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized) from baseline.

- dpam:

  List of survival regressions for model endpoints. This must include
  pre-progression death (PPD), time to progression (TTP) and post
  progression survival calculated under the clock reset model (PPS-CR).

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
prob_pd_stm_cr(0:100, params)
#>   [1] 0.000000000 0.060138735 0.108121914 0.147019738 0.178485112 0.204031115
#>   [7] 0.224898274 0.241978542 0.255814655 0.266761717 0.275114914 0.281129804
#>  [13] 0.285031884 0.287030029 0.287325545 0.286118938 0.283602846 0.279959691
#>  [19] 0.275357940 0.269951702 0.263887342 0.257285072 0.250259636 0.242913869
#>  [25] 0.235336236 0.227605222 0.219789614 0.211948479 0.204132184 0.196383250
#>  [31] 0.188737209 0.181223250 0.173865274 0.166682547 0.159689972 0.152899110
#>  [37] 0.146318278 0.139953175 0.133807059 0.127881545 0.122176413 0.116690118
#>  [43] 0.111419978 0.106362436 0.101513431 0.096867583 0.092419669 0.088163888
#>  [49] 0.084094142 0.080219398 0.076502228 0.072952148 0.069562820 0.066315417
#>  [55] 0.063229318 0.060285215 0.057477267 0.054799591 0.052246535 0.049812641
#>  [61] 0.047492604 0.045281319 0.043173866 0.041165500 0.039251895 0.037428207
#>  [67] 0.035690487 0.034034715 0.032457061 0.030953852 0.029521303 0.028156597
#>  [73] 0.026856198 0.025617453 0.024436951 0.023311968 0.022239959 0.021218396
#>  [79] 0.020244868 0.019317084 0.018434681 0.017591864 0.016788545 0.016022838
#>  [85] 0.015292949 0.014597168 0.013933868 0.013301503 0.012698600 0.012123756
#>  [91] 0.011575636 0.011052971 0.010554551 0.010079221 0.009625886 0.009193506
#>  [97] 0.008781088 0.008387688 0.008012409 0.007654395 0.007312832
# }
```
