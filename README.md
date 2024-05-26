# psm3mkv: A package to evaluate the fit and efficiency of three state oncology cost-effectiveness model structures <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/Merck/psm3mkv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Merck/psm3mkv/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Merck/psm3mkv/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Merck/psm3mkv?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/psm3mkv)](https://cran.r-project.org/package=psm3mkv)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/psm3mkv)](https://cran.r-project.org/package=psm3mkv)
<!-- badges: end -->

The goal of psm3mkv is to evaluate the efficiency and fit of certain
three state model structures to data typical from an oncology clinical
trial, as described in an accompanying article [1]. The package
evaluates the following structures:

- Partitioned Survival Model/analysis (PSM),

- Clock-forward State Transition Model (STM-CF), and

- Clock-reset State Transition Model (STM-CR).

The state transition models differ from each other in that the
transition from progressive disease to death is a function of time from
baseline in the STM-CF and time from progression in the STM-CR [2, 3].

The package requires a patient-level dataset of time to progression
(TTP), progression-free survival (PFS) and overall survival (OS).

Given this, the package enables:

- Fitting a range of models to endpoints relevant each model type:

  - One piece parametric (distributions according to flexsurv).

  - Royston-Parmar splines (1-3 internal knots, hazard/odds/normal
    scales, again as per flexsurv) [4].

  - Two piece parametric (given a time cutoff).

- Selecting the "best fit" survival models for each endpoint (using the
  Akaike Information Criterion, Bayesian Information Criterion or other
  user preference).

- Deriving and presenting likelihoods for the 3 structures so as to
  evaluate fit.

- Presenting the total number of parameters used for each structure, to
  additionally evaluate efficiency.

- Deriving and presenting restricted mean durations by health state and
  for each of the 3 model structures (given a time horizon), to evaluate
  plausibility and structural sensitivity. Additional functionality:

  - Constraining the estimates to ensure that survival is no greater
    than survival in a background lifetable.

  - Applying discounting to obtained discounted restricted means.

  - Bootstrap standard errors can be derived.

- Graphically illustrate observed and fitted membership probabilities,
  to allow visual inspection of fit of the 3 model structures.

Where two piece modeling is used, modelers should be advised to take
care of interpretation and validity in case different cutoff points are
selected for different endpoints.

Additionally, for parametric modeling of STM structures, the model for
survival in the progressive disease state (post progression survival,
PPS) may be a function of an additional arbitrary explanatory variable.
This is intended to enable the exploration of TTP (or some
transformation) as a predictor for PPS.

## Vignettes

The accompanying `vignette("example")` illustrates how the package can
be used for the one-piece parametric and spline modeling.

A second vignette, `vignette("background-mortality")` illustrates how,
after fitting models, estimates of restricted mean durations in health
states can be calculated after constraining for background mortality
from a given life table. Survival is assumed to be no greater than in a
background lifetable.

## Installation

The package requires version R >= 4.1.0 due to use of the native pipe.
Please ensure R is updated first.

### Latest stable release

Install the latest stable release from CRAN:

``` r
install.packages("psm3mkv")
```

### Development version

Install the latest development version from GitHub (this may not be as
stable):

``` r
# install.packages("pak")
pak::pak("Merck/psm3mkv@main")
```

Note that `pak::pak()` does not build the vignettes by default when
installing a package from GitHub, which is ideal because the vignettes
can take a long time to generate. You can conveniently view them on the
package [documentation website](https://merck.github.io/psm3mkv/).

### Additional dependencies

Running the vignettes requires additional dependencies, which are all
either imported by or suggested by *psm3mkv*. Thus you can ensure they
are all installed by specifying `dependencies = TRUE`.

``` r
pak::pak("Merck/psm3mkv@*release", dependencies = TRUE)
```

## Licensing

Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its
affiliates. All rights reserved.

This file is part of the psm3mkv program.

psm3mkv is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

psm3mkv uses third-party R packages which may be distributed under
different licenses.

## References

1. Muston, D. 2024. "Informing Structural Assumptions for Three State Oncology Cost-Effectiveness Models through Model Efficiency and Fit." _Appl Health Econ Health Policy_. DOI: 10.1007/s40258-024-00884-2

2. Jackson, Christopher. 2016. "flexsurv: A Platform for Parametric
   Survival Modeling in R." _Journal of Statistical Software_ 70 (8): 1--33.

3. Woods, Beth S, Eleftherios Sideris, Stephen Palmer, Nick Latimer,
   and Marta Soares. 2020. "Partitioned Survival and State Transition Models
   for Healthcare Decision Making in Oncology: Where Are We Now?"
   _Value in Health_ 23 (12): 1613--1621.

4. Royston, Patrick, and Mahesh KB Parmar. 2002. "Flexible Parametric
   Proportional-Hazards and Proportional-Odds Models for Censored Survival
   Data, with Application to Prognostic Modelling and Estimation of
   Treatment Effects." _Statistics in Medicine_ 21 (15): 2175--2197.
