# Calculate restricted mean durations

Calculates the restricted mean duration, given the form of a parametric
distribution of Royston-Parmar splines

## Usage

``` r
calc_rmd(Tw, type = NA, spec = NA, survobj = NULL)
```

## Arguments

- Tw:

  is the time horizon (weeks) over which the mean should be calculated.

- type:

  is either "par" for regular parametric form (exponential, weibull etc)
  or "spl" for Royston-Parmar splines.

- spec:

  is a list comprising: If type=="par": `dist` is the statistical
  distribution (named per
  [flexsurv::flexsurvreg](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
  and `pars` is a vector of the parameters for that distribution.

  - Exponential distribution (`exp`) requires the rate parameter.

  - Weibull distribution (both `weibullPH` and `weibull` formulations)
    requires the shape and scale parameters.

  - Log-logistic distribution (`llogis`) requires the shape and scale
    parameters.

  - Log-normal distribution (`lnorm`) requires the meanlog and sdlog
    parameters.

  - Gamma and Gompertz distributions (`gamma` and `gompertz`) require
    the shape and rate parameters.

  - Generalized Gamma requires the mu, sigma and Q parameters if using
    the standard parameterization (`gengamma`) or shape, scale and k
    parameters if using the original parameterization (`gengamma.orig`).
    If type=="spl":

  - `gamma` - Vector of parameters describing the baseline spline
    function, as described in
    [flexsurv::flexsurvspline](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
    This may be supplied as a vector with number of elements equal to
    the length of knots, in which case the parameters are common to all
    times. Alternatively a matrix may be supplied, with rows
    corresponding to different times, and columns corresponding to
    knots.

  - `knots` - Vector of locations of knots on the axis of log time,
    supplied in increasing order. Unlike in
    [flexsurv::flexsurvspline](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
    these include the two boundary knots.

  - `scale` - Either "hazard", "odds", or "normal", as described in
    [flexsurv::flexsurvspline](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
    With the default of no knots in addition to the boundaries, this
    model reduces to the Weibull, log-logistic and log-normal
    respectively. The scale must be common to all times.

- survobj:

  is a survival fit object from
  [flexsurv::flexsurvspline](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  or
  [flexsurv::flexsurvreg](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)

## Value

the restricted mean duration, a numeric value.

## Examples

``` r
calc_rmd(Tw=200,
    type="spl",
    spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
    )
#> [1] 199.9296
calc_rmd(Tw=250,
    type="par",
    spec=list(dist="lnorm", pars=c(3,1))
    )
#> [1] 32.45453
```
