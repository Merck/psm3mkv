# Check consistency of PFS definition Check that PFS is defined consistently with TTP and OS in a dataset. This convenience function compares `pfs.durn` with the lower of `ttp.durn` and `os.durn`, and checks that the event field `pfs.flag` is consistent with `ttp.flag` and `os.flag` (is 1 when either `ttp.flag` or `os.flag` is one).

Check consistency of PFS definition Check that PFS is defined
consistently with TTP and OS in a dataset. This convenience function
compares `pfs.durn` with the lower of `ttp.durn` and `os.durn`, and
checks that the event field `pfs.flag` is consistent with `ttp.flag` and
`os.flag` (is 1 when either `ttp.flag` or `os.flag` is one).

## Usage

``` r
check_consistent_pfs(ds)
```

## Arguments

- ds:

  Tibble of complete patient-level dataset

  - `ttp.durn`, `pfs.durn`, and `os.durn` are the durations of TTP (time
    to progression), PFS (progression-free survival), and OS (overall
    survival).

  - `ttp.flag`, `pfs.flag`, and `os.flag`, and `pps.flag` are event flag
    indicators for TTP, PFS, and OS respectively (1=event, 0=censoring).

## Value

List containing:

- `durn`: Logical vector comparing expected and actual PFS durations

- `flag`: Logical vector comparing expected and actual PFS event flags

- `all`: Single logical value of TRUE if all durations and flags match
  as expected, FALSE otherwise

## Examples

``` r
ponc <- create_dummydata("pharmaonc")
check_consistent_pfs(ponc)
#> Warning: Unknown or uninitialised column: `pfs.durn`.
#> $durn
#> logical(0)
#> 
#> $flag
#>   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#>  [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [106] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [121] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [136] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [151] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [166] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [181] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [196] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [211] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [226] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> [241] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#> 
#> $all
#> [1] TRUE
#> 
```
