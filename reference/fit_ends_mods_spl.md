# Fit multiple spline regressions to the multiple required endpoints

Fits multiple survival regressions, according to the distributions
stipulated, to the multiple endpoints required in fitting partitioned
survival analysis, clock forward and clock reset semi-markov models.

## Usage

``` r
fit_ends_mods_spl(
  simdat,
  knot_set = 1:3,
  scale_set = c("hazard", "odds", "normal"),
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

- knot_set:

  is a vector of the numbers of knots to consider, following
  [`flexsurv::flexsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)).

- scale_set:

  is a vector of the spline scales to consider, following
  [`flexsurv::flexsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)).

- expvar:

  Explanatory variable for modeling of PPS

## Value

A list by endpoint, then distribution, each containing two components:

- result: A list of class
  [flexsurv::flexsurvspline](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  containing information about the fitted model.

- error: Any error message returned on fitting the regression (NULL
  indicates no error). Also, the given cuttime.

## See also

Parametric modeling is handled by
[`fit_ends_mods_par()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_par.md)

## Examples

``` r
# \donttest{
# Create dataset in suitable form using bos dataset from the flexsurv package
bosonc <- create_dummydata("flexbosms")
fit_ends_mods_spl(bosonc, expvar=bosonc$ttp.durn)
#> $ttp
#> $ttp[[1]]
#> $ttp[[1]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.7114  -3.1392  -2.2836   0.2183
#> gamma1   0.9650   0.5632   1.3668   0.2050
#> gamma2   0.0110  -0.0211   0.0432   0.0164
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -422.8569, df = 3
#> AIC = 851.7137
#> 
#> 
#> $ttp[[1]]$error
#> NULL
#> 
#> 
#> $ttp[[2]]
#> $ttp[[2]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.8065  -3.2682  -2.3448   0.2356
#> gamma1   0.8694   0.4668   1.2720   0.2054
#> gamma2  -0.0153  -0.0511   0.0205   0.0183
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.7298, df = 3
#> AIC = 849.4597
#> 
#> 
#> $ttp[[2]]$error
#> NULL
#> 
#> 
#> $ttp[[3]]
#> $ttp[[3]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -1.60314  -1.83224  -1.37404   0.11689
#> gamma1   0.38216   0.20690   0.55743   0.08942
#> gamma2  -0.01936  -0.03709  -0.00163   0.00905
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.5655, df = 3
#> AIC = 849.1309
#> 
#> 
#> $ttp[[3]]$error
#> NULL
#> 
#> 
#> $ttp[[4]]
#> $ttp[[4]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -2.73940  -3.27510  -2.20370   0.27332
#> gamma1   0.94983   0.43830   1.46135   0.26098
#> gamma2  -0.00852  -0.14676   0.12973   0.07053
#> gamma3   0.02978  -0.19546   0.25502   0.11492
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -422.8021, df = 4
#> AIC = 853.6043
#> 
#> 
#> $ttp[[4]]$error
#> NULL
#> 
#> 
#> $ttp[[5]]
#> $ttp[[5]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.6568  -3.2137  -2.1000   0.2841
#> gamma1   1.0097   0.4700   1.5495   0.2754
#> gamma2   0.0691  -0.0964   0.2346   0.0844
#> gamma3  -0.1375  -0.4197   0.1448   0.1440
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.3791, df = 4
#> AIC = 850.7583
#> 
#> 
#> $ttp[[5]]$error
#> NULL
#> 
#> 
#> $ttp[[6]]
#> $ttp[[6]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -1.5359  -1.8178  -1.2539   0.1439
#> gamma1   0.4345   0.2068   0.6622   0.1162
#> gamma2   0.0304  -0.0583   0.1190   0.0452
#> gamma3  -0.0796  -0.2359   0.0768   0.0798
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -421.3112, df = 4
#> AIC = 850.6223
#> 
#> 
#> $ttp[[6]]$error
#> NULL
#> 
#> 
#> $ttp[[7]]
#> $ttp[[7]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est        L95%       U95%       se       
#> gamma0  -2.485963  -3.046933  -1.924993   0.286214
#> gamma1   1.197935   0.554813   1.841057   0.328130
#> gamma2   0.221183   0.000829   0.441538   0.112428
#> gamma3  -0.837151  -1.573248  -0.101055   0.375566
#> gamma4   0.887201   0.107971   1.666432   0.397574
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -420.1814, df = 5
#> AIC = 850.3627
#> 
#> 
#> $ttp[[7]]$error
#> NULL
#> 
#> 
#> $ttp[[8]]
#> $ttp[[8]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.4671  -3.0586  -1.8755   0.3018
#> gamma1   1.1985   0.5473   1.8496   0.3322
#> gamma2   0.2335  -0.0139   0.4809   0.1262
#> gamma3  -0.7797  -1.6893   0.1298   0.4641
#> gamma4   0.7355  -0.2786   1.7496   0.5174
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -420.0042, df = 5
#> AIC = 850.0085
#> 
#> 
#> $ttp[[8]]$error
#> NULL
#> 
#> 
#> $ttp[[9]]
#> $ttp[[9]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -1.4492  -1.7515  -1.1470   0.1542
#> gamma1   0.5042   0.2406   0.7677   0.1345
#> gamma2   0.1192  -0.0130   0.2514   0.0674
#> gamma3  -0.4462  -0.9748   0.0824   0.2697
#> gamma4   0.4287  -0.1742   1.0316   0.3076
#> 
#> N = 204,  Events: 103,  Censored: 101
#> Total time at risk: 2358.845
#> Log-likelihood = -420.0067, df = 5
#> AIC = 850.0134
#> 
#> 
#> $ttp[[9]]$error
#> NULL
#> 
#> 
#> 
#> $ppd
#> $ppd[[1]]
#> $ppd[[1]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -5.8254  -7.7259  -3.9250   0.9696
#> gamma1   1.5623   0.4569   2.6676   0.5640
#> gamma2   0.0167  -0.1170   0.1504   0.0682
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.7031, df = 3
#> AIC = 313.4063
#> 
#> 
#> $ppd[[1]]$error
#> NULL
#> 
#> 
#> $ppd[[2]]
#> $ppd[[2]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -5.79819  -7.70574  -3.89064   0.97326
#> gamma1   1.52049   0.39352   2.64746   0.57499
#> gamma2  -0.00446  -0.14762   0.13870   0.07304
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -154.0731, df = 3
#> AIC = 314.1461
#> 
#> 
#> $ppd[[2]]$error
#> NULL
#> 
#> 
#> $ppd[[3]]
#> $ppd[[3]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.7953  -3.4937  -2.0969   0.3563
#> gamma1   0.5645   0.1203   1.0088   0.2267
#> gamma2  -0.0331  -0.0986   0.0324   0.0334
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -154.0678, df = 3
#> AIC = 314.1356
#> 
#> 
#> $ppd[[3]]$error
#> NULL
#> 
#> 
#> $ppd[[4]]
#> $ppd[[4]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -6.260  -8.623  -3.897   1.206
#> gamma1   2.080   0.308   3.852   0.904
#> gamma2   0.404  -0.502   1.309   0.462
#> gamma3  -0.433  -1.434   0.569   0.511
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.3264, df = 4
#> AIC = 314.6528
#> 
#> 
#> $ppd[[4]]$error
#> NULL
#> 
#> 
#> $ppd[[5]]
#> $ppd[[5]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -6.279  -8.668  -3.891   1.219
#> gamma1   2.108   0.301   3.916   0.922
#> gamma2   0.487  -0.491   1.465   0.499
#> gamma3  -0.560  -1.666   0.545   0.564
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.5727, df = 4
#> AIC = 315.1455
#> 
#> 
#> $ppd[[5]]$error
#> NULL
#> 
#> 
#> $ppd[[6]]
#> $ppd[[6]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -2.920  -3.724  -2.115   0.410
#> gamma1   0.757   0.108   1.405   0.331
#> gamma2   0.207  -0.254   0.669   0.236
#> gamma3  -0.287  -0.839   0.266   0.282
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.6044, df = 4
#> AIC = 315.2088
#> 
#> 
#> $ppd[[6]]$error
#> NULL
#> 
#> 
#> $ppd[[7]]
#> $ppd[[7]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -6.3287  -8.9746  -3.6827   1.3500
#> gamma1   2.1917  -0.2409   4.6243   1.2411
#> gamma2   0.2345  -1.2648   1.7337   0.7649
#> gamma3   0.0531  -2.7173   2.8235   1.4135
#> gamma4  -0.3086  -2.3272   1.7101   1.0299
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.3243, df = 5
#> AIC = 316.6486
#> 
#> 
#> $ppd[[7]]$error
#> NULL
#> 
#> 
#> $ppd[[8]]
#> $ppd[[8]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -6.320  -8.969  -3.671   1.352
#> gamma1   2.182  -0.262   4.626   1.247
#> gamma2   0.223  -1.329   1.775   0.792
#> gamma3   0.174  -2.806   3.154   1.521
#> gamma4  -0.487  -2.746   1.772   1.153
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.5704, df = 5
#> AIC = 317.1408
#> 
#> 
#> $ppd[[8]]$error
#> NULL
#> 
#> 
#> $ppd[[9]]
#> $ppd[[9]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -2.89719  -3.73120  -2.06319   0.42552
#> gamma1   0.71931  -0.08609   1.52472   0.41093
#> gamma2   0.00886  -0.65672   0.67443   0.33959
#> gamma3   0.24938  -1.22124   1.71999   0.75033
#> gamma4  -0.38355  -1.60015   0.83305   0.62073
#> 
#> N = 204,  Events: 29,  Censored: 175
#> Total time at risk: 2358.845
#> Log-likelihood = -153.5862, df = 5
#> AIC = 317.1723
#> 
#> 
#> $ppd[[9]]$error
#> NULL
#> 
#> 
#> 
#> $pfs
#> $pfs[[1]]
#> $pfs[[1]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -2.71386  -3.13361  -2.29411   0.21416
#> gamma1   0.99398   0.61132   1.37664   0.19524
#> gamma2   0.00393  -0.01838   0.02623   0.01138
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -512.0123, df = 3
#> AIC = 1030.025
#> 
#> 
#> $pfs[[1]]$error
#> NULL
#> 
#> 
#> $pfs[[2]]
#> $pfs[[2]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est        L95%       U95%       se       
#> gamma0  -2.848058  -3.298634  -2.397482   0.229890
#> gamma1   0.832337   0.449968   1.214706   0.195090
#> gamma2  -0.025453  -0.051358   0.000453   0.013217
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.7148, df = 3
#> AIC = 1029.43
#> 
#> 
#> $pfs[[2]]$error
#> NULL
#> 
#> 
#> $pfs[[3]]
#> $pfs[[3]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est       L95%      U95%      se      
#> gamma0  -1.62401  -1.84681  -1.40122   0.11367
#> gamma1   0.37082   0.20155   0.54008   0.08636
#> gamma2  -0.02243  -0.03539  -0.00947   0.00661
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.4637, df = 3
#> AIC = 1028.927
#> 
#> 
#> $pfs[[3]]$error
#> NULL
#> 
#> 
#> $pfs[[4]]
#> $pfs[[4]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.7744  -3.2640  -2.2848   0.2498
#> gamma1   0.9271   0.4496   1.4046   0.2436
#> gamma2  -0.0307  -0.1621   0.1006   0.0670
#> gamma3   0.0441  -0.1240   0.2122   0.0858
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.8845, df = 4
#> AIC = 1031.769
#> 
#> 
#> $pfs[[4]]$error
#> NULL
#> 
#> 
#> $pfs[[5]]
#> $pfs[[5]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.7070  -3.2250  -2.1890   0.2643
#> gamma1   0.9774   0.4707   1.4841   0.2585
#> gamma2   0.0822  -0.0934   0.2577   0.0896
#> gamma3  -0.1403  -0.3779   0.0973   0.1212
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.2509, df = 4
#> AIC = 1030.502
#> 
#> 
#> $pfs[[5]]$error
#> NULL
#> 
#> 
#> $pfs[[6]]
#> $pfs[[6]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -1.5724  -1.8376  -1.3072   0.1353
#> gamma1   0.4129   0.1989   0.6270   0.1092
#> gamma2   0.0272  -0.0673   0.1217   0.0482
#> gamma3  -0.0637  -0.1942   0.0668   0.0666
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -511.2899, df = 4
#> AIC = 1030.58
#> 
#> 
#> $pfs[[6]]$error
#> NULL
#> 
#> 
#> $pfs[[7]]
#> $pfs[[7]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.5410  -3.0590  -2.0230   0.2643
#> gamma1   1.1468   0.5509   1.7427   0.3040
#> gamma2   0.1899  -0.0307   0.4104   0.1125
#> gamma3  -0.7089  -1.3939  -0.0238   0.3495
#> gamma4   0.6196   0.0378   1.2013   0.2968
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -509.8239, df = 5
#> AIC = 1029.648
#> 
#> 
#> $pfs[[7]]$error
#> NULL
#> 
#> 
#> $pfs[[8]]
#> $pfs[[8]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.5276  -3.0767  -1.9784   0.2802
#> gamma1   1.1501   0.5431   1.7570   0.3097
#> gamma2   0.2371  -0.0218   0.4960   0.1321
#> gamma3  -0.7336  -1.6567   0.1896   0.4710
#> gamma4   0.5475  -0.2959   1.3910   0.4303
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -510.0705, df = 5
#> AIC = 1030.141
#> 
#> 
#> $pfs[[8]]$error
#> NULL
#> 
#> 
#> $pfs[[9]]
#> $pfs[[9]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -1.4793  -1.7636  -1.1950   0.1450
#> gamma1   0.4844   0.2357   0.7330   0.1269
#> gamma2   0.1297  -0.0121   0.2716   0.0724
#> gamma3  -0.4727  -1.0114   0.0660   0.2749
#> gamma4   0.3763  -0.1221   0.8746   0.2543
#> 
#> N = 204,  Events: 132,  Censored: 72
#> Total time at risk: 2358.845
#> Log-likelihood = -509.868, df = 5
#> AIC = 1029.736
#> 
#> 
#> $pfs[[9]]$error
#> NULL
#> 
#> 
#> 
#> $os
#> $os[[1]]
#> $os[[1]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -5.0074  -6.2393  -3.7754   0.6286
#> gamma1   1.3996   0.7437   2.0556   0.3347
#> gamma2  -0.0157  -0.0995   0.0682   0.0428
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -433.2801, df = 3
#> AIC = 872.5603
#> 
#> 
#> $os[[1]]$error
#> NULL
#> 
#> 
#> $os[[2]]
#> $os[[2]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -4.8971  -6.1132  -3.6809   0.6205
#> gamma1   1.2418   0.5761   1.9075   0.3396
#> gamma2  -0.0854  -0.1808   0.0101   0.0487
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -435.5343, df = 3
#> AIC = 877.0686
#> 
#> 
#> $os[[2]]$error
#> NULL
#> 
#> 
#> $os[[3]]
#> $os[[3]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.5243  -3.0238  -2.0248   0.2549
#> gamma1   0.5252   0.2308   0.8195   0.1502
#> gamma2  -0.0771  -0.1250  -0.0292   0.0245
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -434.7993, df = 3
#> AIC = 875.5985
#> 
#> 
#> $os[[3]]$error
#> NULL
#> 
#> 
#> $os[[4]]
#> $os[[4]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -5.321  -6.770  -3.872   0.739
#> gamma1   1.745   0.799   2.692   0.483
#> gamma2   0.299  -0.208   0.807   0.259
#> gamma3  -0.402  -1.050   0.246   0.331
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -432.5721, df = 4
#> AIC = 873.1442
#> 
#> 
#> $os[[4]]$error
#> NULL
#> 
#> 
#> $os[[5]]
#> $os[[5]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -5.369  -6.851  -3.886   0.757
#> gamma1   1.807   0.817   2.797   0.505
#> gamma2   0.614  -0.020   1.249   0.324
#> gamma3  -0.930  -1.784  -0.076   0.436
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -433.4306, df = 4
#> AIC = 874.8611
#> 
#> 
#> $os[[5]]$error
#> NULL
#> 
#> 
#> $os[[6]]
#> $os[[6]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.6769  -3.2520  -2.1018   0.2934
#> gamma1   0.7502   0.3378   1.1625   0.2104
#> gamma2   0.3136  -0.0281   0.6553   0.1743
#> gamma3  -0.5287  -1.0036  -0.0538   0.2423
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -432.8091, df = 4
#> AIC = 873.6182
#> 
#> 
#> $os[[6]]$error
#> NULL
#> 
#> 
#> $os[[7]]
#> $os[[7]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est     L95%    U95%    se    
#> gamma0  -5.426  -6.984  -3.868   0.795
#> gamma1   1.879   0.738   3.019   0.582
#> gamma2   0.318  -0.573   1.209   0.455
#> gamma3  -0.275  -2.322   1.772   1.044
#> gamma4  -0.117  -1.820   1.587   0.869
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -432.4254, df = 5
#> AIC = 874.8508
#> 
#> 
#> $os[[7]]$error
#> NULL
#> 
#> 
#> $os[[8]]
#> $os[[8]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -5.4331  -7.0056  -3.8607   0.8023
#> gamma1   1.8919   0.7269   3.0570   0.5944
#> gamma2   0.3778  -0.6283   1.3840   0.5134
#> gamma3   0.0734  -2.4752   2.6220   1.3003
#> gamma4  -0.8477  -3.1700   1.4745   1.1848
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -433.2146, df = 5
#> AIC = 876.4292
#> 
#> 
#> $os[[8]]$error
#> NULL
#> 
#> 
#> $os[[9]]
#> $os[[9]]$result
#> Call:
#> .f(formula = ..1, k = ..2, scale = ..3)
#> 
#> Estimates: 
#>         est      L95%     U95%     se     
#> gamma0  -2.6970  -3.2959  -2.0981   0.3056
#> gamma1   0.7819   0.3067   1.2570   0.2424
#> gamma2   0.1897  -0.3455   0.7250   0.2731
#> gamma3   0.0191  -1.4401   1.4784   0.7446
#> gamma4  -0.4670  -1.8384   0.9043   0.6997
#> 
#> N = 204,  Events: 97,  Censored: 107
#> Total time at risk: 3479.143
#> Log-likelihood = -432.6368, df = 5
#> AIC = 875.2736
#> 
#> 
#> $os[[9]]$error
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
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[2]]
#> $pps_cf[[2]]$result
#> NULL
#> 
#> $pps_cf[[2]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[3]]
#> $pps_cf[[3]]$result
#> NULL
#> 
#> $pps_cf[[3]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[4]]
#> $pps_cf[[4]]$result
#> NULL
#> 
#> $pps_cf[[4]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[5]]
#> $pps_cf[[5]]$result
#> NULL
#> 
#> $pps_cf[[5]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[6]]
#> $pps_cf[[6]]$result
#> NULL
#> 
#> $pps_cf[[6]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[7]]
#> $pps_cf[[7]]$result
#> NULL
#> 
#> $pps_cf[[7]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[8]]
#> $pps_cf[[8]]$result
#> NULL
#> 
#> $pps_cf[[8]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cf[[9]]
#> $pps_cf[[9]]$result
#> NULL
#> 
#> $pps_cf[[9]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> 
#> $pps_cr
#> $pps_cr[[1]]
#> $pps_cr[[1]]$result
#> NULL
#> 
#> $pps_cr[[1]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[2]]
#> $pps_cr[[2]]$result
#> NULL
#> 
#> $pps_cr[[2]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[3]]
#> $pps_cr[[3]]$result
#> NULL
#> 
#> $pps_cr[[3]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[4]]
#> $pps_cr[[4]]$result
#> NULL
#> 
#> $pps_cr[[4]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[5]]
#> $pps_cr[[5]]$result
#> NULL
#> 
#> $pps_cr[[5]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[6]]
#> $pps_cr[[6]]$result
#> NULL
#> 
#> $pps_cr[[6]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[7]]
#> $pps_cr[[7]]$result
#> NULL
#> 
#> $pps_cr[[7]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[8]]
#> $pps_cr[[8]]$result
#> NULL
#> 
#> $pps_cr[[8]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> $pps_cr[[9]]
#> $pps_cr[[9]]$result
#> NULL
#> 
#> $pps_cr[[9]]$error
#> <simpleError in model.frame.default(formula = survival::Surv(time = durn1, time2 = durn2,     event = evflag) ~ {    {        expvar    }}, data = <environment>): variable lengths differ (found for '{     {         expvar     } }')>
#> 
#> 
#> 
# }
```
