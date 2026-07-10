# Calculate likelihoods for three three-state model structures

Calculate likelihood values and other summary output for the following
three state models structures: partitioned survival, clock forward state
transition, and clock reset state transition. The function requires
appropriately formatted patient-level data, a set of fitted survival
regressions, and the time cut-off (if two-piece modeling is used).

## Usage

``` r
calc_likes(ptdata, dpam, cuttime = 0)
```

## Arguments

- ptdata:

  Dataset of patient level data. Must be a tibble with columns named:

  - `ptid`: patient identifier

  - `pfs.durn`: duration of PFS from baseline

  - `pfs.flag`: event flag for PFS (=1 if progression or death occurred,
    0 for censoring)

  - `os.durn`: duration of OS from baseline

  - `os.flag`: event flag for OS (=1 if death occurred, 0 for censoring)

  - `ttp.durn`: duration of TTP from baseline (usually should be equal
    to pfs.durn)

  - `ttp.flag`: event flag for TTP (=1 if progression occurred, 0 for
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

- cuttime:

  Time cutoff - this is nonzero for two-piece models.

## Value

A list of three tibbles: `all` is a tibble of results for all patients:

- `methname`: the model structure or method.

- `npar`: is the number of parameters used by that method.

- `npts_1` to `npts_4` are the number of patients experiencing outcomes
  1-4 respectively (see below), and `npts_tot` the total.

- `ll_1` to `ll_4` are the log-likelihood values for patients
  experiencing outcomes 1-4 respectively (see below), and `ll_tot` the
  total. `valid` is a tibble of the same design as `all` but only in
  patients with valid likelihoods for all 4 methods `sum` is a tibble in
  respect of patients with valid likelihoods for all 4 methods
  providing:

- `npts`: number of patients contributing results for this method.

- `npar`: number of parameters used by that method.

- `ll`: total log-likelihood

- `AIC`: Akaike Information Criterion value for this model

- `BIC`: Bayesian Information Criterion value for this model

The four outcomes are as follows:

- \(1\) refers to patients who remain alive and progression-free during
  the follow-up;

- \(2\) refers to patients who die without prior progression during the
  follow-up;

- \(3\) refers to patients who progress and then remain alive for the
  remaining follow-up, and

- \(4\) refers to patients who progress and die within the follow-up.

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
calc_likes(bosonc, dpam=params)
#> $all
#> # A tibble: 4 × 12
#>   methname     npar npts_1 npts_2 npts_3 npts_4 npts_tot  ll_1  ll_2  ll_3  ll_4
#>   <chr>       <dbl>  <int>  <int>  <int>  <int>    <int> <dbl> <dbl> <dbl> <dbl>
#> 1 psm_simple      7     72     29     35     68      204 -64.8 -152. -151.   NA 
#> 2 psm_complex     9     72     29     35     68      204 -64.8 -147. -150. -473.
#> 3 stm_cf          9     72     29     35     68      204 -65.2 -147. -148. -473.
#> 4 stm_cr         10     72     29     35     68      204 -65.2 -147. -148. -470.
#> # ℹ 1 more variable: ll_tot <dbl>
#> 
#> $valid
#> # A tibble: 4 × 12
#>   methname     npar npts_1 npts_2 npts_3 npts_4 npts_tot  ll_1  ll_2  ll_3  ll_4
#>   <chr>       <dbl>  <int>  <int>  <int>  <int>    <int> <dbl> <dbl> <dbl> <dbl>
#> 1 psm_simple      7     72     29     35     67      203 -64.8 -152. -151. -469.
#> 2 psm_complex     9     72     29     35     67      203 -64.8 -147. -150. -468.
#> 3 stm_cf          9     72     29     35     67      203 -65.2 -147. -148. -468.
#> 4 stm_cr         10     72     29     35     67      203 -65.2 -147. -148. -464.
#> # ℹ 1 more variable: ll_tot <dbl>
#> 
#> $sum
#> # A tibble: 4 × 8
#>   methname     npts  npar    ll   aic   bic rank_aic rank_bic
#>   <chr>       <int> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 psm_simple    203     7 -836. 1685. 1709.        4        4
#> 2 psm_complex   203     9 -830. 1679. 1708.        3        3
#> 3 stm_cf        203     9 -828. 1673. 1703.        2        2
#> 4 stm_cr        203    10 -824. 1668. 1701.        1        1
#> 
# }
```
