# Derive PPS survival function under a PSM

Derive the post-progression survival (PPS) function under the simple or
complex PSM formulation.

## Usage

``` r
calc_surv_psmpps(totime, fromtime = 0, ptdata, dpam, psmtype = "simple")
```

## Arguments

- totime:

  Vector of times to which the survival function is calculated

- fromtime:

  Vector of times from which the survival function is calculated

- ptdata:

  Patient-level dataset

- dpam:

  List of fitted survival models for each endpoint

- psmtype:

  Either "simple" or "complex" PSM formulation

## Value

Vector of PPS survival function values

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
calc_surv_psmpps(totime=1:10,
  fromtime=rep(1,10),
  ptdata=bosonc,
  dpam=params,
  psmtype="simple")
#>  [1] 1.0000000 0.9995962 0.9833657 0.9541989 0.9190001 0.8807816 0.8411186
#>  [8] 0.8009357 0.7607707 0.7209705
# }
```
