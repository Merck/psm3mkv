# Mortality adjustments

## Introduction

This is the second in a series of vignettes illustrating methods for
evaluating the fit and efficiency of three state oncology
cost-effectiveness model structures, as described in an accompanying
journal article (Muston 2024). The package is heavily dependent on
flexsurv (Jackson 2016).

After fitting models, as described in
[`vignette("example")`](https://merck.github.io/psm3mkv/articles/example.md),
estimates of Restricted Mean Durations (RMDs) in health states can be
calculated after constraining for background mortality from a given life
table.

### Notation

We denote functions of survival, hazard and cumulative hazard as
$`S(t), h(t)`$ and $`H(t)`$. There are ‘unadjusted’ values for relevant
endpoints (e.g. PFS) dependent on the model structure (PSM, STM-CF or
STM-CR). We assume that there is a ‘general background mortality’,
against which any extrapolations and derivations of RMD should be
constrained, denoted with subscript ‘gen’. Consequently there are
‘adjusted’ values, denoted with superscript ‘adj’.

### Background mortality survival and hazard functions

Our lifetable gives us decreasing $`l_x`$ values at points between
$`t=0, ..., t_{max}`$, with $`l_{t_{max}}=0`$. Let us define the
background mortality survival function as follows.

``` math
S_{gen}(t) = 1 - R \left(1 - \frac{l_{x_0 + t}}{l_{x_0}} \right)
```

where $`x_0`$ is the mean age of the study population at baseline, $`R`$
is a [Standardized Mortality
Ratio](https://en.wikipedia.org/wiki/Standardized_mortality_ratio), and
$`l_x`$ is an appropriate lifetable.

We can also derive the average background mortality hazard for $`t`$ in
the range $`[t_1, t_2)`$.

``` math
h_{gen}(t) = \frac{\log[S_{gen}(t_1)] - \log[S_{gen}(t_2)]}{t_2 - t_1}
```

## What does constraining for background mortality mean?

There are a number of approaches applied in literature and HTA
submissions, so it could mean several things!

### Approach 1 - Constraining only the survival function

Under this approach, the background mortality is assumed to limit the
survival function as follows. There is no direct constraint to the
hazard.

``` math
S^{adj}(t) = \min \left[ S(t), S_{gen}(t) \right]
```

This is a simple approach, and has been used at least one NICE
technology appraisal (National Institute for Health and Care Excellence
2022):

> “Patient longevity is always the lesser of values generated from the
> disease-specific survival curve (after adjustment for treatment and
> functional status) and the survival curve for the general population
> according to age and sex.”

However, this approach does not ensure that the hazard is at least as
great as the hazard of background mortality.

### Approach 2 - Constraining the hazard function

A more rigorous adjustment would be that the hazard of background
mortality acts as a constraint to the unadjusted hazard (Sweeting et al.
2023).

``` math
h^{adj}(t) = \max \left[ h(t), h_{gen}(t) \right]
```
``` math
\implies S^{adj}(t) = \exp \left[ - \int_0^t h^{adj}(u) du \right]
```

In either case, we must assume that the original dataset was NOT subject
to background mortality. Although unlikely to be true, this is a common,
pragmatic assumption in cost-effectiveness models as long as background
mortality is relatively insignificant during trial follow-up.

### Approach 3 - Modeling excess hazard in the original dataset

The above approaches work by adjusting extrapolations made from models
fitted to data assumed to be not subject to background mortality. If -
as is likely - the population in the dataset were in fact subject to
background mortality, then it would be better to model that dataset
using excess hazard methods rather than to seek to make adjustments to
extrapolations after the fact. Further discussion of survival
extrapolation incorporating general population mortality is provided by
Sweeting et al. (2023).

## Application to 3-state State Transition Models (STMs)

As discussed previously, state transition models with three states may
be clock forward (‘STM-CF’) or ‘clock reset’ (STM-CR). Such structures
involve modeling each transition directly, which may be performed by
fitting survival distributions to the TTP, PPD and PPS endpoints.

In general, we wish to constrain the pre- and post-progression mortality
(PPD and PPS) are at least as great as background mortality.

\$\$ h^{adj}\_{PPD}(t) = \max \left\[h\_{gen}(t), h\_{PPD}(t) \right\]
\\ h^{adj}\_{PPS}(t) = \max \left\[h\_{gen}(t), h\_{PPS}(t) \right\] \\
\$\$

Formulas for the mean time in PF and OS are shown in the accompanying
article in respect of hazard and survival functions for TTP, PPD and PPS
endpoints, where PPS is a function of time from baseline in a Clock
Forward (CF) and time from progression in a Clock Reset (CR) model.
Adjusting for background mortality requires using the adjusted rather
than unadjusted hazards for PPD and PPS in formulae.

This is straightforward enough for the PPD endpoint in the STM-CF and
STM-CR structures, and for the PPS endpoint in the STM-CF structure.
However matters are more complex for the PPS endpoint for the STM-CF
since its hazard and survival are functions of two times (time from
progression and time from baseline) rather than just one (time from
baseline). This package does not provide integral solutions relying on
continuous time for RMDs, and must instead rely on discretization.

Spreadsheet-based economic models almost always rely on discretizing
time into timesteps rather than follow the integral formulas described
above for continuous time. As long as timesteps are reasonably short,
restricted mean duration results should remain reasonably accurate.
Derided as a “kludge”, half-cycle corrections are nevertheless
recommended (Naimark et al. 2013). The package provides discretized
calculations.

## Application to Partitioned Survival Models (PSM)

The RMD in the PF state and time alive are integrals over the time
horizon, $`T`$, of the survival functions of the endpoints PFS and OS.
Endpoints PPD and PPS are not explicitly defined.

However as discussed in the accompanying manuscript, TTP - and therefore
the PPD and PPS endpoints which directly follow once PFS and OS are
defined - are often in practice implicitly defined in PSM economic model
in order to calculate payoffs that depend on progression events. Two
methods are possible:

- Simple (default): the hazard of TTP relative to PFS is proportionate
  to the progression events observed over the follow-up period relative
  to the total number of observed PFS events.

- Complex: the hazard of TTP is directly modeled by fitting survival
  distributions as with other endpoints.

The adjustment to RMD estimates can then proceed as with the STMs, with
hazards and survival functions of TTP, PPD and PPS that reflect the PSM
structure (in which PFS and OS endpoints are first estimated), a chosen
model of progression (simple or complex), and then finally the
adjustment to PPD and PPS for background mortality.

## Illustration of calculations

### Set-up

First we load the packages we need - all of which are suggested for or
imported to *psm3mkv*. If you haven’t installed *psm3mkv* yet, please
see the [installation
instructions](https://github.com/Merck/psm3mkv?tab=readme-ov-file#installation)
to install it with its dependencies. Of these, we will making further
use of Wickham et al. (2023), Riffe (2015) and Müller and Wickham
(2023).

``` r

library("dplyr")
library("psm3mkv")
library("tibble")
```

### Creating the life table

In order to apply constraints to background mortality, we need some
background mortality data, in the form of a lifetable. We can take this
from the Human Mortality Database using the *HMDHFDplus* package
mentioned above (with thanks to Robert Hettle for the recommendation).
The lifetable will need to start from an assumed age at baseline.
Lifetables are constructed by time in years.

``` r

# HMD HFD package
library("HMDHFDplus")

# Assumed population age at baseline (time=0)
baseage <- 51.0

# Mortality data - England & Wales, Female, 2019
mort <- readHMDweb(
  CNTRY = "GBRTENW",
  item = "fltper_1x1",
  username = "",
  password = ""
) |>
  filter(Year == 2019) |>
  select(Age, lx, dx) |>
  mutate(Timey = ceiling(Age - baseage)) |>
  filter(Timey >= 0)

# The table needs to end with lx=0
mort <- add_row(mort, Age = 111, lx = 0, Timey = 60)
```

You will see that the above code cannot run without a login for the
[Human Mortality Database](https://www.mortality.org). Alternatively, we
could just make up a mortality table.

``` r

mort <- tibble(
  Timey = 0:30,
  lx = 10000 * exp(-0.03 * Timey^1.1),
  dx = lx - lead(lx),
  qx = dx / lx
)
```

However we do it, once we have the mortality data, we can apply the SMR
($`R`$) and derive a lifetable from time zero as required.

``` r

# Assumed Standardized Mortality Ratio
SMR <- 2

# Recalculate the lifetable with the SMR applied
mort$adjlx <- mort$lx
for (i in 2:length(mort$lx)) {
  mort$adjlx[i] <- mort$adjlx[i - 1] * (1 - SMR * mort$qx[i - 1])
}

# Ensure lx>=0
mort$adjlx[mort$adjlx < 0] <- 0

# Create and view lifetable
ltable <- tibble(lttime = mort$Timey, lx = mort$adjlx)
head(ltable)
#> # A tibble: 6 × 2
#>   lttime     lx
#>    <int>  <dbl>
#> 1      0 10000 
#> 2      1  9409.
#> 3      2  8774.
#> 4      3  8151.
#> 5      4  7553.
#> 6      5  6985.
```

### Obtaining and fitting distributions to patient-level data

Next we get some patient-level data and fit the PSM and STMs. So far, we
are closely following
[`vignette("example")`](https://merck.github.io/psm3mkv/articles/example.md).
Time is assumed to be recorded in the patient-level data in weeks, with
approximately 52.18 weeks per year.

``` r

# Get some data
bosonc <- create_dummydata("flexbosms")

# We'll make the durations a lot longer
# so that they will be definitely constrained by the lifetable
bosonc <- bosonc |>
  mutate(
    pfs.durn = 20 * pfs.durn,
    os.durn = 20 * os.durn,
    ttp.durn = 20 * ttp.durn
  )

# Fit parametric models to each endpoint
allfits <- fit_ends_mods_par(bosonc)

# Bring together preferred fits for each endpoint
params <- list(
  ppd = find_bestfit(allfits$ppd, "aic")$fit,
  ttp = find_bestfit(allfits$ttp, "aic")$fit,
  pfs = find_bestfit(allfits$pfs, "aic")$fit,
  os = find_bestfit(allfits$os, "aic")$fit,
  pps_cf = find_bestfit(allfits$pps_cf, "aic")$fit,
  pps_cr = find_bestfit(allfits$pps_cr, "aic")$fit
)
```

### Making the projections

We are skipping over any internal or external validation for the
purposes of this vignette, and jump instead straight to making the
survival projections. We will assume a 10 year time horizon.

First we derive a projection using integral formulae, without lifetable
constraints.

``` r

# Set time horizon
thoz <- 10

# Run the calculations
proj1 <- calc_allrmds(bosonc, dpam = params, Ty = thoz)

# Present the results
res1 <- proj1$results |> mutate(method = "int", lxadj = "no", psmmeth = "simple")
res1 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)
#> # A tibble: 3 × 7
#>          pf        pd        os model  method lxadj psmmeth
#>   <num:.1!> <num:.1!> <num:.1!> <chr>  <chr>  <chr> <chr>  
#> 1     274.4     118.3     392.7 PSM    int    no    simple 
#> 2     270.2     119.7     389.9 STM-CF int    no    simple 
#> 3     270.2     123.8     393.9 STM-CR int    no    simple
```

Next we derive a projection using a discretization approximation, still
without any lifetable constraints.

``` r

# Run the calculation
proj2 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, rmdmethod = "disc")

# Present the results
res2 <- proj2$results |>
  mutate(method = "disc", lxadj = "no", psmmeth = "simple")
res2 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)
#> # A tibble: 3 × 7
#>          pf        pd        os model  method lxadj psmmeth
#>   <num:.1!> <num:.1!> <num:.1!> <chr>  <chr>  <chr> <chr>  
#> 1     274.3     118.3     392.6 PSM    disc   no    simple 
#> 2     270.1     119.7     389.8 STM-CF disc   no    simple 
#> 3     270.1     123.7     393.8 STM-CR disc   no    simple
```

Discretization has proven fairly accurate. The STM-CR estimate of mean
time alive has reduced by just 0.1 weeks for example, from 393.9 to
393.8 weeks.

Finally, we use a discretization approximation and apply the lifetable
constraints.

``` r

# Run the calculations
proj3 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, rmdmethod = "disc", lifetable = ltable)

# Present the results
res3 <- proj3$results |>
  mutate(method = "disc", lxadj = "yes", psmmeth = "simple")
res3 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)
#> # A tibble: 3 × 7
#>          pf        pd        os model  method lxadj psmmeth
#>   <num:.1!> <num:.1!> <num:.1!> <chr>  <chr>  <chr> <chr>  
#> 1     235.4     105.5     340.9 PSM    disc   yes   simple 
#> 2     230.1     107.6     337.7 STM-CF disc   yes   simple 
#> 3     230.1     111.6     341.6 STM-CR disc   yes   simple

# Proportion of time alive spent in PF
proppf3 <- res3$pf / res3$os
```

There is now some difference between the STM-CR compared to the STM-CF
and PSM structures. All structures estimate fairly similar RMD estimates
of being alive (range: 337.7 to 341.6 weeks). The estimate of mean time
alive in the STM-CR structure has reduced by 52.3 weeks from 393.9 to
341.6 weeks.

The proportion of that time in the PF rather than PD state varies
somewhat (PSM: 69%, STM-CF: 68.1%, STM-CR: 67.3%). Given that the hazard
of mortality would typically be lower before rather than after
progression, the pre-progression mortality is likely to be more affected
by any mortality constraint that post-progression mortality. The STM-CR
model of post-progression survival depends on time from progression and
is therefore rather different from the STM-CF and PSM structures, which
model mortality based on time from baseline.

We can check whether using a ‘complex’ type of PSM structure rather than
the default ‘simple’ affects RMD estimates from the PSM structure.

``` r

# Run the calculations
proj4 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, psmtype = "complex", rmdmethod = "disc", lifetable = ltable)

# Present the results
res4 <- proj4$results |>
  mutate(method = "disc", lxadj = "yes", psmmeth = "complex")
res4 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)
#> # A tibble: 3 × 7
#>          pf        pd        os model  method lxadj psmmeth
#>   <num:.1!> <num:.1!> <num:.1!> <chr>  <chr>  <chr> <chr>  
#> 1     231.0     103.7     334.7 PSM    disc   yes   complex
#> 2     230.1     107.6     337.7 STM-CF disc   yes   complex
#> 3     230.1     111.6     341.6 STM-CR disc   yes   complex

# Proportion of time alive spent in PF
proppf4 <- res4$pf / res4$os
```

All structures estimate fairly similar RMD estimates of being alive
(range: 334.7 to 341.6 weeks). The proportion of that time in the PF
rather than PD state varied a little by model structure (PSM: 69%,
STM-CF: 68.1%, STM-CR: 67.3%).

### Comparing the results

A summary of the STM-CF estimates of mean time in PF state and mean time
alive is given in the table below.

| Model | Mean time in PF (weeks) | Change (weeks) | Mean time alive (weeks) | Change (weeks) |
|:---|:---|:---|:---|:---|
| No lifetable constraint, integral methods | 270.2 | \- | 389.9 | \- |
| No lifetable constraint, discretized calculation | 270.1 | -0.1 | 389.8 | -0.1 |
| With lifetable constraint, discretized calculation | 230.1 | -40.1 | 337.7 | -52.2 |

In this fictional case, the application of a lifetable constraint had
fairly sizeable effects on the results - although this was by design of
this demonstration.

A summary of the STM-CR estimates of mean time in PF state and mean time
alive is given in the table below.

| Model | Mean time in PF (weeks) | Change (weeks) | Mean time alive (weeks) | Change (weeks) |
|:---|:---|:---|:---|:---|
| No lifetable constraint, integral methods | 270.2 | \- | 393.9 | \- |
| No lifetable constraint, discretized calculation | 270.1 | -0.1 | 393.8 | -0.1 |
| With lifetable constraint, discretized calculation | 230.1 | -40.1 | 341.6 | -52.3 |

The application of a lifetable constraint had similar and fairly
sizeable effects on the results here also.

A summary of the PSM estimates of mean time in PF state and mean time
alive is given in the table below.

| Model | Mean time in PF (weeks) | Change (weeks) | Mean time alive (weeks) | Change (weeks) |
|:---|:---|:---|:---|:---|
| No lifetable constraint, integral methods, simple PSM | 274.4 | \- | 392.7 | \- |
| No lifetable constraint, discretized calculation, simple PSM | 274.3 | -0.1 | 392.6 | -0.1 |
| With lifetable constraint, discretized calculation, simple PSM | 235.4 | -39 | 340.9 | -51.8 |
| With lifetable constraint, discretized calculation, complex PSM | 231 | -43.4 | 334.7 | -58.1 |

The change in the RMDs for the ‘complex’ PSM structure compared to the
‘simple’ PSM structure were -4.3, -1.9 and -6.2 weeks respectively for
time in PF state, PD state and time alive.

### Possible extensions

Not applied here, the package also allows application of discounting
through the `discrate` optional call to `calc_allrnds()`. This work
could also be easily extended through calculating RMDs or Life Years
that are Quality Adjusted (QALYs).

## References

Jackson, Christopher. 2016. “flexsurv: A Platform for Parametric
Survival Modeling in R.” *Journal of Statistical Software* 70 (8): 1–33.

Müller, Kirill, and Hadley Wickham. 2023. *tibble: Simple Data Frames*.
<https://CRAN.R-project.org/package=tibble>.

Muston, Dominic. 2024. “Informing Structural Assumptions for Three State
Oncology Cost-Effectiveness Models Through Model Efficiency and Fit.”
*Applied Health Economics and Health Policy*, ahead of print.
<https://doi.org/10.1007/s40258-024-00884-2>.

Naimark, David MJ, Nader N Kabboul, and Murray D Krahn. 2013. “The
Half-Cycle Correction Revisited: Redemption of a Kludge.” *Medical
Decision Making* 33 (7): 961–70.

National Institute for Health and Care Excellence. 2022.
*Avalglucosidase Alfa for Treating Pompe Disease*. Technology Appraisal
Guidance No. TA821. National Institute for Health; Care Excellence.
<https://www.nice.org.uk/guidance/ta821/documents/committee-papers>.

Riffe, Tim. 2015. *Reading Human Fertility Database and Human Mortality
Database Data into R*. TR-2015-004. MPIDR.
<https://doi.org/10.4054/MPIDR-TR-2015-004>.

Sweeting, Michael J, Mark J Rutherford, Dan Jackson, et al. 2023.
“Survival Extrapolation Incorporating General Population Mortality Using
Excess Hazard and Cure Models: A Tutorial.” *Medical Decision Making* 43
(6): 737–48.

Wickham, Hadley, Romain François, Lionel Henry, Kirill Müller, and Davis
Vaughan. 2023. *dplyr: A Grammar of Data Manipulation*.
<https://CRAN.R-project.org/package=dplyr>.
