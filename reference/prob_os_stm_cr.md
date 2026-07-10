# Calculate probability of being alive under the state transition clock reset model

Calculates membership probability of being alive at a given time
(vectorized). This probability is from the state transition clock reset
model, according to the given statistical distributions and parameters.

## Usage

``` r
prob_os_stm_cr(time, dpam, starting = c(1, 0, 0))
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
prob_os_stm_cr(0:100, params)
#>   [1] 1.00000000 0.99695476 0.98950782 0.97746752 0.96161310 0.94299857
#>   [7] 0.92254723 0.90090563 0.87846191 0.85545262 0.83203515 0.80832010
#>  [13] 0.78438687 0.76030439 0.73614308 0.71198223 0.68790283 0.66398489
#>  [19] 0.64030304 0.61692553 0.59392015 0.57133517 0.54921944 0.52761537
#>  [25] 0.50655596 0.48606814 0.46617256 0.44688418 0.42821302 0.41016450
#>  [31] 0.39274008 0.37593750 0.35975163 0.34417487 0.32919723 0.31480714
#>  [37] 0.30099142 0.28773575 0.27502473 0.26284252 0.25117260 0.23999818
#>  [43] 0.22930225 0.21906785 0.20927831 0.19991639 0.19096573 0.18241001
#>  [49] 0.17423324 0.16643504 0.15896917 0.15183686 0.14502391 0.13850405
#>  [55] 0.13228962 0.12635463 0.12068691 0.11527456 0.11010623 0.10517103
#>  [61] 0.10045851 0.09595866 0.09166190 0.08755905 0.08364158 0.07990060
#>  [67] 0.07632835 0.07291716 0.06965971 0.06654902 0.06357815 0.06074126
#>  [73] 0.05803193 0.05544476 0.05297372 0.05061357 0.04835939 0.04620636
#>  [79] 0.04414988 0.04218558 0.04031110 0.03851873 0.03680657 0.03517099
#>  [85] 0.03360852 0.03211587 0.03068988 0.02932754 0.02802600 0.02678250
#>  [91] 0.02559444 0.02445933 0.02337478 0.02233852 0.02134838 0.02040230
#>  [97] 0.01949830 0.01863449 0.01780908 0.01702035 0.01626666
# }
```
