# VLOOKUP function

Function to lookup values according to an index. Aims to behave
similarly to VLOOKUP in Microsoft Excel, however several lookups can be
made at once (`indexval` can be a vector) and interpolation is available
where lookups are inexact (choice of 4 methods).

## Usage

``` r
vlookup(indexval, indexvec, valvec, method = "geom")
```

## Arguments

- indexval:

  The index value to be looked-up (may be a vector of multiple values)

- indexvec:

  The vector of indices to look-up within

- valvec:

  The vector of values corresponding to the vector of indices

- method:

  Method may be `floor`, `ceiling`, `arith` or `geom` (default).

## Value

Numeric value or vector, depending on the lookup/interpolation method
chosen:

- `floor`: Floor (minimum) value, where interpolation is required
  between measured values

- `ceiling`: Ceiling (maximum) value, where interpolation is required
  between measured values

- `arith`: Arithmetic mean, where interpolation is required between
  measured values

- `geom`: Geometric mean, where interpolation is required between
  measured values

## See also

[`HMDHFDplus::readHMDweb()`](https://rdrr.io/pkg/HMDHFDplus/man/readHMDweb.html)
can be used to obtain lifetables from the Human Mortality Database

## Examples

``` r
# Suppose we have survival probabilities at times 0 to 20
times <- 0:20
survival <- 1-times*0.04
# We would like to look-up the survival probability at time 7
vlookup(7, times, survival)
#> [1] 0.72
# In this case, the floor, ceiling, arith and geom values are identical
# because survival time 7 is known, and no interpolation is necessary
vlookup(c(7, 7.5), times, survival)
#> [1] 0.7200000 0.6997142
# The second row of the returned tibble reveal different estimates of the survival at time 7.5.
# The values vary according to the interpolation method between
# observed survival values at times 7 and 8.
```
