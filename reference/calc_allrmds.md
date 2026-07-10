# Calculate restricted mean durations for each health state and all three models

Calculate restricted mean durations for each health state (progression
free and progressed disease) for all three models (partitioned survival,
clock forward state transition model, clock reset state transition
model).

## Usage

``` r
calc_allrmds(
  ptdata,
  inclset = 0,
  dpam,
  psmtype = "simple",
  cuttime = 0,
  Ty = 10,
  lifetable = NA,
  discrate = 0,
  rmdmethod = "int",
  timestep = 1,
  boot = FALSE
)
```

## Arguments

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

- inclset:

  Vector to indicate which patients to include in analysis

- dpam:

  List of statistical fits to each endpoint required in PSM, STM-CF and
  STM-CR models.

- psmtype:

  Either "simple" or "complex" PSM formulation

- cuttime:

  Time cutoff - this is nonzero for two-piece models.

- Ty:

  Time duration over which to calculate. Assumes input is in years, and
  patient-level data is recorded in weeks.

- lifetable:

  Optional, a life table. Columns must include `lttime` (time in years,
  or 52.18 times shorter than the time index elsewhere, starting from
  zero) and `lx`

- discrate:

  Discount rate (% per year)

- rmdmethod:

  can be "int" (default for full integral calculations) or "disc" for
  approximate discretized calculations

- timestep:

  required if method=="int", default being 1

- boot:

  logical flag to indicate whether abbreviated output is required
  (default = FALSE), for example for bootstrapping

## Value

List of detailed numeric results

- cutadj indicates the survival function and area under the curves for
  PFS and OS up to the cutpoint

- results provides results of the restricted means calculations, by
  model and state.

## Examples

``` r
# \donttest{
# Create dataset and fit survival models (splines)
bosonc <- create_dummydata("flexbosms")
fits <- fit_ends_mods_par(bosonc)
# Pick out best distribution according to min AIC
params <- list(
  ppd = find_bestfit(fits$ppd, "aic")$fit,
  ttp = find_bestfit(fits$ttp, "aic")$fit,
  pfs = find_bestfit(fits$pfs, "aic")$fit,
  os = find_bestfit(fits$os, "aic")$fit,
  pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
  pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
)
# RMD using default "int" method, no lifetable constraint
calc_allrmds(bosonc, dpam=params)
#> $cutadj
#> $cutadj$pfarea
#> [1] 0
#> 
#> $cutadj$pfsurv
#> [1] 1
#> 
#> $cutadj$osarea
#> [1] 0
#> 
#> $cutadj$ossurv
#> [1] 1
#> 
#> 
#> $results
#> # A tibble: 3 × 4
#>      pf    pd    os model 
#>   <dbl> <dbl> <dbl> <chr> 
#> 1  17.9  9.82  27.7 PSM   
#> 2  20.0 10.7   30.7 STM-CF
#> 3  20.0 10.1   30.0 STM-CR
#> 
# RMD using discretized ("disc") method, no lifetable constraint
calc_allrmds(bosonc, dpam=params, rmdmethod="disc", timestep=1, boot=TRUE)
#> [1] 17.87470 19.96263 19.96263 11.98304 10.73811 10.07684 29.85774 30.70075
#> [9] 30.03948
# }
```
