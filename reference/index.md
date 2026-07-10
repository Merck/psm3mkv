# Package index

## All functions

- [`calc_allrmds()`](https://merck.github.io/psm3mkv/reference/calc_allrmds.md)
  : Calculate restricted mean durations for each health state and all
  three models

- [`calc_haz_psm()`](https://merck.github.io/psm3mkv/reference/calc_haz_psm.md)
  : Derive pre and post-progression hazards of death under PSM

- [`calc_likes()`](https://merck.github.io/psm3mkv/reference/calc_likes.md)
  : Calculate likelihoods for three three-state model structures

- [`calc_rmd()`](https://merck.github.io/psm3mkv/reference/calc_rmd.md)
  : Calculate restricted mean durations

- [`calc_surv_psmpps()`](https://merck.github.io/psm3mkv/reference/calc_surv_psmpps.md)
  : Derive PPS survival function under a PSM

- [`check_consistent_pfs()`](https://merck.github.io/psm3mkv/reference/check_consistent_pfs.md)
  :

  Check consistency of PFS definition Check that PFS is defined
  consistently with TTP and OS in a dataset. This convenience function
  compares `pfs.durn` with the lower of `ttp.durn` and `os.durn`, and
  checks that the event field `pfs.flag` is consistent with `ttp.flag`
  and `os.flag` (is 1 when either `ttp.flag` or `os.flag` is one).

- [`compare_psm_likes()`](https://merck.github.io/psm3mkv/reference/compare_psm_likes.md)
  : Compare likelihoods of PSMs

- [`constrain_survprob()`](https://merck.github.io/psm3mkv/reference/constrain_survprob.md)
  : Constrain survival probabilities according to hazards in a lifetable
  Recalculated constrained survival probabilities (by week) as the lower
  of the original unadjusted survival probability and the survival
  implied by the given lifetable (assumed indexed as years).

- [`create_dummydata()`](https://merck.github.io/psm3mkv/reference/create_dummydata.md)
  : Create dummy dataset for illustration

- [`create_extrafields()`](https://merck.github.io/psm3mkv/reference/create_extrafields.md)
  : Create the additional time-to-event endpoints, adjusting for
  cutpoint

- [`find_bestfit()`](https://merck.github.io/psm3mkv/reference/find_bestfit.md)
  : Find the "best" survival regression from a list of model fits

- [`fit_ends_mods_par()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_par.md)
  : Fit multiple parametric survival regressions to the multiple
  required endpoints

- [`fit_ends_mods_spl()`](https://merck.github.io/psm3mkv/reference/fit_ends_mods_spl.md)
  : Fit multiple spline regressions to the multiple required endpoints

- [`graph_psm_hazards()`](https://merck.github.io/psm3mkv/reference/graph_psm_hazards.md)
  : Graph the PSM hazard functions

- [`graph_psm_survs()`](https://merck.github.io/psm3mkv/reference/graph_psm_survs.md)
  : Graph the PSM survival functions

- [`graph_survs()`](https://merck.github.io/psm3mkv/reference/graph_survs.md)
  : Graph the observed and fitted state membership probabilities

- [`prob_os_psm()`](https://merck.github.io/psm3mkv/reference/prob_os_psm.md)
  : Calculate probability of being alive in a partitioned survival model

- [`prob_os_stm_cf()`](https://merck.github.io/psm3mkv/reference/prob_os_stm_cf.md)
  : Calculate probability of being alive under the state transition
  clock forward model

- [`prob_os_stm_cr()`](https://merck.github.io/psm3mkv/reference/prob_os_stm_cr.md)
  : Calculate probability of being alive under the state transition
  clock reset model

- [`prob_pd_psm()`](https://merck.github.io/psm3mkv/reference/prob_pd_psm.md)
  : Calculate membership probability of progressed disease state in a
  partitioned survival model

- [`prob_pd_stm_cf()`](https://merck.github.io/psm3mkv/reference/prob_pd_stm_cf.md)
  : Calculate probability of having progressed disease under the state
  transition clock forward model

- [`prob_pd_stm_cr()`](https://merck.github.io/psm3mkv/reference/prob_pd_stm_cr.md)
  : Calculate probability of having progressed disease under the state
  transition clock reset model

- [`prob_pf_psm()`](https://merck.github.io/psm3mkv/reference/prob_pf_psm.md)
  : Calculate probability of being progression free in partitioned
  survival model

- [`prob_pf_stm()`](https://merck.github.io/psm3mkv/reference/prob_pf_stm.md)
  : Calculate probability of being progression free in either state
  transition model (clock forward or clock reset)

- [`prob_pps_cf()`](https://merck.github.io/psm3mkv/reference/prob_pps_cf.md)
  : Calculate probability of post progression survival under the state
  transition clock forward model

- [`prob_pps_cr()`](https://merck.github.io/psm3mkv/reference/prob_pps_cr.md)
  : Calculate probability of post progression survival under the state
  transition clock reset model

- [`vlookup()`](https://merck.github.io/psm3mkv/reference/vlookup.md) :
  VLOOKUP function
