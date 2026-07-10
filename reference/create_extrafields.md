# Create the additional time-to-event endpoints, adjusting for cutpoint

Create the additional time-to-event endpoints, adjusting for cutpoint

## Usage

``` r
create_extrafields(ds, cuttime = 0)
```

## Arguments

- ds:

  Patient-level dataset

- cuttime:

  Time cutpoint

## Value

Tibble of complete patient-level dataset, adjusted for cutpoint
`ttp.durn`, `pfs.durn`, `ppd.durn` and `os.durn` are the durations of
TTP (time to progression), PFS (progression-free survival), PPD
(pre-progression death) and OS (overall survival) respectively beyond
the cutpoint. `pps.durn` is the duration of survival beyond progression,
irrespective of the cutpoint. `pps.odurn` is the difference between
`ttp.durn` and `os.durn` (which may be different to `pps.durn`).
`ttp.flag`, `pfs.flag`, `ppd.flag`, `os.flag`, and `pps.flag` are event
flag indicators for TTP, PFS, PPD, OS and PPS respectively (1=event,
0=censoring).

## Examples

``` r
bosonc <- create_dummydata("flexbosms")
create_extrafields(bosonc, cuttime=10)
#> # A tibble: 204 × 16
#>     ptid pfs.odurn pfs.flag os.odurn os.flag ttp.odurn ttp.flag ttp.durn
#>    <int>     <dbl>    <dbl>    <dbl>   <dbl>     <dbl>    <dbl>    <dbl>
#>  1     1    18.7          1    42.9        1    18.7          1    8.67 
#>  2     2    12.0          1    23.3        1    12.0          1    2.01 
#>  3     3     0.452        1     8.81       1     0.452        1    0    
#>  4     4     9.07         1    52.7        1     9.07         1    0    
#>  5     5    47.7          0    47.7        0    47.7          0   37.7  
#>  6     6     3.26         1    13.1        1     3.26         1    0    
#>  7     7     0.738        1     2.74       1     0.738        1    0    
#>  8     8    35.5          1    35.5        1    35.5          0   25.5  
#>  9     9    10.4          1    37.1        1    10.4          1    0.369
#> 10    10    10.5          1    32.7        1    10.5          1    0.464
#> # ℹ 194 more rows
#> # ℹ 8 more variables: pfs.durn <dbl>, os.durn <dbl>, tzero <dbl>,
#> #   pps.odurn <dbl>, pps.durn <dbl>, pps.flag <dbl>, ppd.durn <dbl>,
#> #   ppd.flag <dbl>
```
