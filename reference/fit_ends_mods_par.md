# Fit multiple parametric survival regressions to the multiple required endpoints

Fits multiple parametric survival regressions, according to the
distributions stipulated, to the multiple endpoints required in fitting
partitioned survival analysis, clock forward and clock reset semi-markov
models.

## Usage

``` r
fit_ends_mods_par(
  simdat,
  cuttime = 0,
  ppd.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  ttp.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  pfs.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  os.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  pps_cf.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  pps_cr.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
  expvar = NA
)
```

## Arguments

- simdat:

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

- cuttime:

  Cut-off time for a two-piece model, equals zero for one-piece models.

- ppd.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Pre-Progression Death (PPD).

- ttp.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Time To Progression (TTP).

- pfs.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Progression-Free Survival (PFS).

- os.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Overall Survival (OS).

- pps_cf.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Post Progression Survival, where time is from baseline
  (clock forward).

- pps_cr.dist:

  Vector of distributions (named per
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  to be fitted to Post Progression Survival, where time is from
  progression (clock reset).

- expvar:

  Explanatory variable for modeling of PPS

## Value

A list by endpoint, then distribution, each containing two components:

- result: A list of class *flexsurvreg* containing information about the
  fitted model.

- error: Any error message returned on fitting the regression (NULL
  indicates no error).

## See also

Spline modeling is handled by
[`fit_ends_mods_spl()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_spl.md)

## Examples

``` r
bosonc <- create_dummydata("flexbosms")
fit_ends_mods_par(bosonc, expvar=bosonc$ttp.durn)
#> $ttp
#> $ttp[[1]]
#> $ttp[[1]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>       est     L95%    U95%    se    
#> rate  0.0437  0.0360  0.0530  0.0043
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -425.5134, df = 1
#> AIC = 853.0269
#> 
#> 
#> $ttp[[1]]$error
#> NULL
#> 
#> 
#> $ttp[[2]]
#> $ttp[[2]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape  0.8379  0.7110  0.9875  0.0702
#> scale  0.0680  0.0448  0.1033  0.0145
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -423.0927, df = 2
#> AIC = 850.1853
#> 
#> 
#> $ttp[[2]]$error
#> NULL
#> 
#> 
#> $ttp[[3]]
#> $ttp[[3]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape   1.0300   0.8726   1.2158   0.0872
#> scale  15.4175  11.8479  20.0625   2.0716
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -422.068, df = 2
#> AIC = 848.136
#> 
#> 
#> $ttp[[3]]$error
#> NULL
#> 
#> 
#> $ttp[[4]]
#> $ttp[[4]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>          est    L95%   U95%   se   
#> meanlog  2.762  2.462  3.062  0.153
#> sdlog    1.777  1.535  2.057  0.133
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -423.7432, df = 2
#> AIC = 851.4864
#> 
#> 
#> $ttp[[4]]$error
#> NULL
#> 
#> 
#> $ttp[[5]]
#> $ttp[[5]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  0.80871  0.65294  1.00164  0.08828
#> rate   0.03166  0.02120  0.04730  0.00649
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -423.5205, df = 2
#> AIC = 851.0411
#> 
#> 
#> $ttp[[5]]$error
#> NULL
#> 
#> 
#> $ttp[[6]]
#> $ttp[[6]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  -0.0376  -0.0652  -0.0101   0.0141
#> rate    0.0603   0.0455   0.0800   0.0087
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.3123, df = 2
#> AIC = 846.6246
#> 
#> 
#> $ttp[[6]]$error
#> NULL
#> 
#> 
#> 
#> $ppd
#> $ppd[[1]]
#> $ppd[[1]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>       est      L95%     U95%     se     
#> rate  0.01229  0.00854  0.01769  0.00228
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -156.5603, df = 1
#> AIC = 315.1206
#> 
#> 
#> $ppd[[1]]$error
#> NULL
#> 
#> 
#> $ppd[[2]]
#> $ppd[[2]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  1.43500  1.09372  1.88278  0.19884
#> scale  0.00352  0.00105  0.01183  0.00218
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.7336, df = 2
#> AIC = 311.4673
#> 
#> 
#> $ppd[[2]]$error
#> NULL
#> 
#> 
#> $ppd[[3]]
#> $ppd[[3]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape   1.553   1.174   2.055   0.222
#> scale  43.036  29.365  63.072   8.393
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -154.0749, df = 2
#> AIC = 312.1498
#> 
#> 
#> $ppd[[3]]$error
#> NULL
#> 
#> 
#> $ppd[[4]]
#> $ppd[[4]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>          est    L95%   U95%   se   
#> meanlog  3.923  3.450  4.397  0.242
#> sdlog    1.289  0.993  1.674  0.172
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -154.541, df = 2
#> AIC = 313.0821
#> 
#> 
#> $ppd[[4]]$error
#> NULL
#> 
#> 
#> $ppd[[5]]
#> $ppd[[5]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape  1.5585  1.0961  2.2159  0.2799
#> rate   0.0303  0.0152  0.0605  0.0107
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.7877, df = 2
#> AIC = 311.5754
#> 
#> 
#> $ppd[[5]]$error
#> NULL
#> 
#> 
#> $ppd[[6]]
#> $ppd[[6]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  0.03933  0.00680  0.07185  0.01659
#> rate   0.00777  0.00436  0.01385  0.00229
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -154.1277, df = 2
#> AIC = 312.2554
#> 
#> 
#> $ppd[[6]]$error
#> NULL
#> 
#> 
#> 
#> $pfs
#> $pfs[[1]]
#> $pfs[[1]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>       est      L95%     U95%     se     
#> rate  0.05596  0.04718  0.06637  0.00487
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -512.5726, df = 1
#> AIC = 1027.145
#> 
#> 
#> $pfs[[1]]$error
#> NULL
#> 
#> 
#> $pfs[[2]]
#> $pfs[[2]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape  0.9313  0.8080  1.0733  0.0675
#> scale  0.0676  0.0453  0.1009  0.0138
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -512.0729, df = 2
#> AIC = 1028.146
#> 
#> 
#> $pfs[[2]]$error
#> NULL
#> 
#> 
#> $pfs[[3]]
#> $pfs[[3]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape   1.1921   1.0310   1.3784   0.0883
#> scale  11.4917   9.3244  14.1628   1.2254
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -513.4304, df = 2
#> AIC = 1030.861
#> 
#> 
#> $pfs[[3]]$error
#> NULL
#> 
#> 
#> $pfs[[4]]
#> $pfs[[4]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>          est    L95%   U95%   se   
#> meanlog  2.412  2.179  2.645  0.119
#> sdlog    1.535  1.353  1.742  0.099
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -516.8304, df = 2
#> AIC = 1037.661
#> 
#> 
#> $pfs[[4]]$error
#> NULL
#> 
#> 
#> $pfs[[5]]
#> $pfs[[5]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  0.91518  0.75087  1.11544  0.09240
#> rate   0.04985  0.03633  0.06840  0.00805
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -512.1786, df = 2
#> AIC = 1028.357
#> 
#> 
#> $pfs[[5]]$error
#> NULL
#> 
#> 
#> $pfs[[6]]
#> $pfs[[6]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est       L95%      U95%      se      
#> shape  -0.01338  -0.03435   0.00759   0.01070
#> rate    0.06342   0.04926   0.08165   0.00817
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.7442, df = 2
#> AIC = 1027.488
#> 
#> 
#> $pfs[[6]]$error
#> NULL
#> 
#> 
#> 
#> $os
#> $os[[1]]
#> $os[[1]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>       est      L95%     U95%     se     
#> rate  0.02788  0.02285  0.03402  0.00283
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -444.2435, df = 1
#> AIC = 890.4871
#> 
#> 
#> $os[[1]]$error
#> NULL
#> 
#> 
#> $os[[2]]
#> $os[[2]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  1.51432  1.29287  1.77369  0.12215
#> scale  0.00559  0.00254  0.01233  0.00226
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -433.3483, df = 2
#> AIC = 870.6965
#> 
#> 
#> $os[[2]]$error
#> NULL
#> 
#> 
#> $os[[3]]
#> $os[[3]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape   1.814   1.541   2.136   0.151
#> scale  23.832  20.361  27.895   1.914
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -436.9493, df = 2
#> AIC = 877.8987
#> 
#> 
#> $os[[3]]$error
#> NULL
#> 
#> 
#> $os[[4]]
#> $os[[4]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>          est     L95%    U95%    se    
#> meanlog  3.1627  2.9845  3.3409  0.0909
#> sdlog    1.0166  0.8806  1.1734  0.0744
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -439.3941, df = 2
#> AIC = 882.7881
#> 
#> 
#> $os[[4]]$error
#> NULL
#> 
#> 
#> $os[[5]]
#> $os[[5]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est     L95%    U95%    se    
#> shape  1.7620  1.3981  2.2207  0.2080
#> rate   0.0605  0.0435  0.0842  0.0102
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -434.1308, df = 2
#> AIC = 872.2617
#> 
#> 
#> $os[[5]]$error
#> NULL
#> 
#> 
#> $os[[6]]
#> $os[[6]]$result
#> Call:
#> .f(formula = ..1, dist = ..2)
#> 
#> Estimates: 
#>        est      L95%     U95%     se     
#> shape  0.04179  0.02473  0.05885  0.00871
#> rate   0.01532  0.01079  0.02174  0.00274
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -433.9344, df = 2
#> AIC = 871.8687
#> 
#> 
#> $os[[6]]$error
#> NULL
#> 
#> 
#> 
#> $pps_cf
#> $pps_cf[[1]]
#> $pps_cf[[1]]$result
#> NULL
#> 
#> $pps_cf[[1]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[2]]
#> $pps_cf[[2]]$result
#> NULL
#> 
#> $pps_cf[[2]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[3]]
#> $pps_cf[[3]]$result
#> NULL
#> 
#> $pps_cf[[3]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[4]]
#> $pps_cf[[4]]$result
#> NULL
#> 
#> $pps_cf[[4]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[5]]
#> $pps_cf[[5]]$result
#> NULL
#> 
#> $pps_cf[[5]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[6]]
#> $pps_cf[[6]]$result
#> NULL
#> 
#> $pps_cf[[6]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> 
#> $pps_cr
#> $pps_cr[[1]]
#> $pps_cr[[1]]$result
#> NULL
#> 
#> $pps_cr[[1]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[2]]
#> $pps_cr[[2]]$result
#> NULL
#> 
#> $pps_cr[[2]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[3]]
#> $pps_cr[[3]]$result
#> NULL
#> 
#> $pps_cr[[3]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[4]]
#> $pps_cr[[4]]$result
#> NULL
#> 
#> $pps_cr[[4]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[5]]
#> $pps_cr[[5]]$result
#> NULL
#> 
#> $pps_cr[[5]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[6]]
#> $pps_cr[[6]]$result
#> NULL
#> 
#> $pps_cr[[6]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> 
```
