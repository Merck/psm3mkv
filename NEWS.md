# psm3mkv (development version)

# psm3mkv 0.3.1

- Submission to CRAN, including changes requested by CRAN

# psm3mkv 0.3.0

- First submission to CRAN, not accepted

# psm3mkv 0.2.2 (4 May 2024)

Several minor changes to ready the package for CRAN.

- Reduced exported functions to only those necessary.
- Reviewed and updated documentation and vignettes.
- Reworked `calc_surv()` and `calc_haz()` functions and calls to make more efficient use of `flexsurv` objects.
- Switched logo to corporate teal.
- Fix all `R CMD check` notes.

# psm3mkv 0.2.1 (14 Apr 2024)

Constraining calculations of restricted mean durations and the accompanying `vignette("background-mortality")` has been reworked. The calculations using integral/continuous methods was not reliable. Instead, `calc_allrmds` now has a `rmdmethod="disc"` option to allow for discretized calculations for a given timestep (defaulting at one week). There is a new collection of functions in `discrmd.R` to provide for this, as well as `constrain_survprob()` function, which constrains a vector of survival estimates at given times such that the underlying hazard is at least as great as in an accompanying lifetable.

# psm3mkv 0.2.0 (26 Jan 2024)

This version provides additional functionality to the calculation of restricted mean durations in `calc_allrmds`. These estimates may now be constrained by a lifetable (see `calc_ltsurv`) and discounting may now be applied. A vignette describing how to use this functionality is provided: `vignette("background-mortality")`.

There are also some changes to _pkgdown.yml to avoid the theming changes imposed in the 0.6.0 release of bslib 0.6.0, with thanks to @nanxstats.

# psm3mkv 0.1.1 (5 Jan 2024)

## Added PSM analysis functions (experimental)

Merged some experimental functions into the main branch following the version 0.1 release package. These functions provide analyses of the constraints on mortality hazards and therefore survival implied by a PSM. They are: `calc_haz_psm()`, `calc_surv_psmpps()`, `pickout_psmhaz`, `graph_psm_hazards()`, and `graph_psm_survs()`.

# psm3mkv 0.1.0 (1 Jan 2024)

This is the initial release of the package, rather belatedly. The code dates to October 2023.
