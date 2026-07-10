# Calculate probability of being alive in a partitioned survival model

Calculates membership probability of being alive at a particular time
(vectorized), given either state transition model (clock forward or
clock reset) with given statistical distributions and parameters. This
is the sum of membership probabilities in the progression free and
progressed disease states.

## Usage

``` r
prob_os_psm(time, dpam, starting = c(1, 0, 0))
```

## Arguments

- time:

  Time (numeric and vectorized)

- dpam:

  List of survival regressions for model endpoints. This must include
  overall survival (OS).

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
prob_os_psm(0:100, params)
#>   [1] 1.000000000 0.993333811 0.982460429 0.969074887 0.953734324 0.936776991
#>   [7] 0.918449615 0.898950119 0.878446448 0.857086233 0.835002270 0.812315818
#>  [13] 0.789138714 0.765574757 0.741720655 0.717667076 0.693502704 0.669312345
#>  [19] 0.645174506 0.621161570 0.597339994 0.573770507 0.550508301 0.527603233
#>  [25] 0.505100023 0.483038462 0.461453625 0.440376084 0.419832131 0.399843999
#>  [31] 0.380430085 0.361605177 0.343380677 0.325764826 0.308762924 0.292377547
#>  [37] 0.276608761 0.261454330 0.246909917 0.232969284 0.219624475 0.206866001
#>  [43] 0.194683013 0.183063466 0.171994277 0.161461472 0.151450328 0.141945500
#>  [49] 0.132931146 0.124391040 0.116308675 0.108667364 0.101450325 0.094640765
#>  [55] 0.088221917 0.082177105 0.076489851 0.071143939 0.066123454 0.061412820
#>  [61] 0.056996834 0.052860693 0.048990017 0.045370867 0.041989763 0.038833690
#>  [67] 0.035890110 0.033146965 0.030592679 0.028216158 0.026006789 0.023954430
#>  [73] 0.022049410 0.020282517 0.018644991 0.017128511 0.015725185 0.014427539
#>  [79] 0.013228500 0.012121389 0.011099902 0.010158095 0.009290377 0.008491485
#>  [85] 0.007756480 0.007080727 0.006459880 0.005889871 0.005366897 0.004887403
#>  [91] 0.004448071 0.004045807 0.003677730 0.003341156 0.003033593 0.002752723
#>  [97] 0.002496397 0.002262621 0.002049549 0.001855470 0.001678803
# }
```
