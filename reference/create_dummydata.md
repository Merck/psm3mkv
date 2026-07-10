# Create dummy dataset for illustration

Create dummy dataset to illustrate
[psm3mkv](https://merck.github.io/psm3mkv/reference/psm3mkv-package.md)

## Usage

``` r
create_dummydata(dsname)
```

## Arguments

- dsname:

  Dataset name, as follows:

  - `flexbosms` provides a dataset based on
    [`flexsurv::bosms3()`](http://chjackson.github.io/flexsurv-dev/reference/bos.md).
    This contains all the fields necessary for
    [psm3mkv](https://merck.github.io/psm3mkv/reference/psm3mkv-package.md).
    Durations have been converted from months in the original dataset to
    weeks.

  - `pharmaonc` provides a dataset based on
    [pharmaverseadam::adsl](https://pharmaverse.github.io/pharmaverseadam/reference/adsl.html)
    and
    [pharmaverseadam::adrs_onco](https://pharmaverse.github.io/pharmaverseadam/reference/adrs_onco.html)
    to demonstrate how this package can be used with ADaM ADTTE
    datasets.

  - `survcan` provides a dataset based on
    [`survival::cancer()`](https://rdrr.io/pkg/survival/man/lung.html).
    This contains the necessary ID and overall survival fields only.
    Durations have been converted from days in the original dataset to
    weeks. You will additionally need to supply PFS and TTP data (fields
    pfs.durn, pfs.flag, ttp.durn and ttp.flag) to use
    [psm3mkv](https://merck.github.io/psm3mkv/reference/psm3mkv-package.md).

## Value

Tibble dataset, for use with
[psm3mkv](https://merck.github.io/psm3mkv/reference/psm3mkv-package.md)
functions

## Examples

``` r
create_dummydata("survcan") |> head()
#>   ptid   os.durn os.flag
#> 1    1  43.71429       1
#> 2    2  65.00000       1
#> 3    3 144.28571       0
#> 4    4  30.00000       1
#> 5    5 126.14286       1
#> 6    6 146.00000       0
create_dummydata("flexbosms") |> head()
#> # A tibble: 6 × 7
#>    ptid pfs.durn pfs.flag os.durn os.flag ttp.durn ttp.flag
#>   <int>    <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
#> 1     1   18.7          1   42.9        1   18.7          1
#> 2     2   12.0          1   23.3        1   12.0          1
#> 3     3    0.452        1    8.81       1    0.452        1
#> 4     4    9.07         1   52.7        1    9.07         1
#> 5     5   47.7          0   47.7        0   47.7          0
#> 6     6    3.26         1   13.1        1    3.26         1
create_dummydata("pharmaonc") |> head()
#> # A tibble: 6 × 7
#>   ptid        os.durn pfs.dur ttp.durn os.flag pfs.flag ttp.flag
#>   <chr>         <dbl>   <dbl>    <dbl>   <dbl>    <dbl>    <dbl>
#> 1 01-701-1015   25.9        9        9       0        0        0
#> 2 01-701-1023    4          0        0       0        0        0
#> 3 01-701-1028   25.6        6        6       0        1        1
#> 4 01-701-1033    3.86       0        0       0        0        0
#> 5 01-701-1034   26          6        6       0        0        0
#> 6 01-701-1047    7.71       0        0       0        0        0
```
