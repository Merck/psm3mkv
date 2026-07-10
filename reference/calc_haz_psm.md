# Derive pre and post-progression hazards of death under PSM

Derive the hazards of death pre- and post-progression under either
simple or complex PSM formulations.

## Usage

``` r
calc_haz_psm(timevar, ptdata, dpam, psmtype)
```

## Arguments

- timevar:

  Vector of times at which to calculate the hazards

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

  Survival data for all other endpoints (time to progression,
  pre-progression death, post-progression survival) are derived from PFS
  and OS.

- dpam:

  List of survival regressions for each endpoint:

  - pre-progression death (PPD)

  - time to progression (TTP)

  - progression-free survival (PFS)

  - overall survival (OS)

  - post-progression survival clock forward (PPS-CF) and

  - post-progression survival clock reset (PPS-CR).

- psmtype:

  Either "simple" or "complex" PSM formulation

## Value

List of pre, the pre-progression hazard, and post, the post-progression
hazard

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
calc_haz_psm(0:10, ptdata=bosonc, dpam=params, psmtype="simple")
#> $adj
#> $adj$ttp
#>  [1] 0.00000000 0.04742078 0.04722152 0.04707233 0.04682576 0.04650723
#>  [7] 0.04614387 0.04575442 0.04533051 0.04479742 0.04418546
#> 
#> $adj$ppd
#>  [1] 0.000000000 0.009903692 0.013295380 0.013253374 0.013183952 0.013094269
#>  [7] 0.012991963 0.012882314 0.012762960 0.012612865 0.012440568
#> 
#> $adj$pfs
#>  [1] 0.00000000 0.05732447 0.06051690 0.06032570 0.06000971 0.05960150
#>  [7] 0.05913583 0.05863674 0.05809348 0.05741028 0.05662603
#> 
#> $adj$os
#>  [1] 0.000000000 0.009300758 0.012253823 0.014430772 0.016195615 0.017678549
#>  [7] 0.018943502 0.020027650 0.020955175 0.021743175 0.022404577
#> 
#> $adj$pps
#>  [1] 0.000000000 0.000000000 0.005094092 0.024838109 0.034454269 0.040318847
#>  [7] 0.044430336 0.047603965 0.050244788 0.052614795 0.054815057
#> 
#> 
#> $unadj
#> $unadj$ttp
#>  [1] 0.00000000 0.04742078 0.04722152 0.04707233 0.04682576 0.04650723
#>  [7] 0.04614387 0.04575442 0.04533051 0.04479742 0.04418546
#> 
#> $unadj$ppd
#>  [1] 0.00000000 0.01335148 0.01329538 0.01325337 0.01318395 0.01309427
#>  [7] 0.01299196 0.01288231 0.01276296 0.01261286 0.01244057
#> 
#> $unadj$pfs
#>  [1] 0.00000000 0.06077226 0.06051690 0.06032570 0.06000971 0.05960150
#>  [7] 0.05913583 0.05863674 0.05809348 0.05741028 0.05662603
#> 
#> $unadj$os
#>  [1] 0.000000000 0.009363175 0.012472587 0.014891286 0.016981265 0.018871673
#>  [7] 0.020625521 0.022278934 0.023854812 0.025368713 0.026831756
#> 
#> $unadj$pps
#>  [1]         NaN 0.000000000 0.005094092 0.024838109 0.034454269 0.040318847
#>  [7] 0.044430336 0.047603965 0.050244788 0.052614795 0.054815057
#> 
#> 
calc_haz_psm(0:10, ptdata=bosonc, dpam=params, psmtype="complex")
#> $adj
#> $adj$ttp
#>  [1] 0.00000000 0.05758869 0.05337416 0.05075486 0.04878425 0.04717919
#>  [7] 0.04581438 0.04459666 0.04337732 0.04216807 0.04098853
#> 
#> $adj$ppd
#>  [1] 0.000000000 0.003183577 0.007142746 0.009570846 0.011225462 0.012422312
#>  [7] 0.013321453 0.014040074 0.014716157 0.015242215 0.015637501
#> 
#> $adj$pfs
#>  [1] 0.00000000 0.06077226 0.06051690 0.06032570 0.06000971 0.05960150
#>  [7] 0.05913583 0.05863674 0.05809348 0.05741028 0.05662603
#> 
#> $adj$os
#>  [1] 0.000000000 0.009300758 0.012253823 0.014430772 0.016195615 0.017678549
#>  [7] 0.018943502 0.020027650 0.020955175 0.021743175 0.022404577
#> 
#> $adj$pps
#>  [1] 0.00000000 0.11640992 0.06026855 0.04720160 0.04346609 0.04281332
#>  [7] 0.04340284 0.04448366 0.04559770 0.04699859 0.04859870
#> 
#> 
#> $unadj
#> $unadj$ttp
#>  [1] 0.00000000 0.05758869 0.05337416 0.05075486 0.04878425 0.04717919
#>  [7] 0.04581438 0.04459666 0.04337732 0.04216807 0.04098853
#> 
#> $unadj$ppd
#>  [1] 0.000000000 0.003183577 0.007142746 0.009570846 0.011225462 0.012422312
#>  [7] 0.013321453 0.014040074 0.014716157 0.015242215 0.015637501
#> 
#> $unadj$pfs
#>  [1] 0.00000000 0.06077226 0.06051690 0.06032570 0.06000971 0.05960150
#>  [7] 0.05913583 0.05863674 0.05809348 0.05741028 0.05662603
#> 
#> $unadj$os
#>  [1] 0.000000000 0.009363175 0.012472587 0.014891286 0.016981265 0.018871673
#>  [7] 0.020625521 0.022278934 0.023854812 0.025368713 0.026831756
#> 
#> $unadj$pps
#>  [1]        NaN 0.11640992 0.06026855 0.04720160 0.04346609 0.04281332
#>  [7] 0.04340284 0.04448366 0.04559770 0.04699859 0.04859870
#> 
#> 
# }
```
