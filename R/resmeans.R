#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the psm3mkv program.
#
#  psm3mkv is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==================================================================
# Functions that calculate restricted mean durations (RMDs)
# resmeans.R
# =========================================================

# Calculation of RMDs by model and state
#   rmd_pf_stm: PF state, according to STM-CR or STM-CF models
#   rmd_pd_stm_cr: PD state, according to STM-CR model
#   rmd_pd_stm_cf: PD state, according to STM-CF model
#   rmd_pf_psm: PF state, according to PSM model
#   rmd_os_psm: OS states, according to PSM model
# Functions that bring RMD calculation results together
#   calc_allrmds: presents means in all states by all models
#   calc_allrmds_boot: version of calc_means presenting limited output for bootstrapping analyses

#' Restricted mean duration in progression-free for state transition models
#' @description Calculates the mean duration in the progression-free state for both the state transition clock forward and clock reset models. Requires a carefully formatted list of fitted survival regressions for the necessary endpoints, and the time duration to calculate over.
#' @param dpam List of survival regressions for model endpoints. These must include time to progression (TTP) and pre-progression death (PPD).
#' @param Ty Time duration over which to calculate. Assumes input is in years, and patient-level data is recorded in weeks.
#' @param starting Vector of membership probabilities at time zero.
#' @param lifetable Optional. The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param discrate Discount rate (%) per year
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
#' @seealso Used safely as [prmd_pf_stm] by [calc_allrmds]
#' @export
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' rmd_pf_stm(dpam=params)
rmd_pf_stm <- function(dpam, Ty=10, starting=c(1, 0, 0), lifetable=NA, discrate=0) {
  # Declare local variables
  Tw <- ttp.ts <- ppd.ts <- NULL
  # Time horizon in weeks
  Tw <- convert_yrs2wks(Ty)
  # Normalize starting vector
  starting <- starting/sum(starting)
  # Pull out type and spec for TTP and PPD
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ppd.ts <- convert_fit2spec(dpam$ppd)
  # RMD PFS is the integral of S_PFS; S_PFS = S_TTP x S_PPD
  # subject to a maximum of the lifetable survival
  integrand <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x))
    sttp <- calc_surv(x, ttp.ts$type, ttp.ts$spec)
    sppd <- calc_surv(x, ppd.ts$type, ppd.ts$spec)
    if (!is.data.frame(lifetable)) {
      vn*sttp*sppd
    } else {
      vn*sttp*pmin(sppd, calc_ltsurv(convert_wks2yrs(x), lifetable))
    }
  }
  int <- stats::integrate(integrand, 0, Tw)
  return(starting[1]*int$value)
}

#' Safely calculate restricted mean duration in progression-free for state transition models
#' @description Calculates the mean duration in the progression-free state for both the state transition clock forward and clock reset models. Requires a carefully formatted list of fitted survival regressions for the necessary endpoints, and the time duration to calculate over. Wrapper with 'possibly' of [rmd_pf_stm]. This function is called by [calc_allrmds].
#' @param ... Pass-through to [rmd_pf_stm]
#' @include basics.R
#' @return Numeric value in same time unit as patient-level data (weeks).
prmd_pf_stm <- purrr::possibly(rmd_pf_stm, otherwise=NA_real_)

#' Restricted mean duration in progressed disease state for clock reset state transition model
#' @description Calculates the mean duration in the progressed disease state for the clock reset state transition model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over.
#' @inheritParams rmd_pf_stm
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @seealso [rmd_pd_stm_cr]
#' @include basics.R
#' @export
#' @seealso Used safely as [prmd_pd_stm_cr] by [calc_allrmds]
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' rmd_pd_stm_cr(dpam=params)
rmd_pd_stm_cr <- function(dpam, Ty=10, starting=c(1, 0, 0), lifetable=NA, discrate=0) {
  # Declare error if starting[2]>0 and lifetable!=NA
  if ((is.data.frame(lifetable))*starting[2]!=0) stop("Cannot apply lifetable constraint when PD state not initially empty")
  # Declare local variables
  Tw <- ttp.ts <- ppd.ts <- pps.ts <- NULL
  S <- int_pf <- int_pd <- soj <- NULL
  # Bound to aid integration in weeks
  Tw <- convert_yrs2wks(Ty)
  # Normalize starting vector
  starting <- starting/sum(starting)
  # Pull out type and spec for TTP, PPD and PPS_CR
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ppd.ts <- convert_fit2spec(dpam$ppd)
  pps.ts <- convert_fit2spec(dpam$pps_cr)
  # Integrand from PF = S_TTP(x1) * S_PPD(x1) * h_TTP(x1) * S_PPS(x2-x1)
  # = S_PPD x f_TTP x S_PPS
  # subject to a maximum of the lifetable survival
  integrand_pf <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x[2]))
    sppd <- calc_surv(x[1], ppd.ts$type, ppd.ts$spec)
    fttp <- calc_dens(x[1], ttp.ts$type, ttp.ts$spec)
    spps <- calc_surv(x[2]-x[1], pps.ts$type, pps.ts$spec)
    if (is.data.frame(lifetable)) {
      gens1 <- calc_ltsurv(convert_wks2yrs(x[1]), lifetable)
      gens2 <- calc_ltsurv(convert_wks2yrs(x[2]), lifetable)
      sppd <- pmin(sppd, gens1)
      spps <- pmin(spps, gens2/gens1)
    }
    vn*sppd*fttp*spps
  }
  S <- cbind(c(0,0),c(0, Tw),c(Tw, Tw))
  int_pf <- SimplicialCubature::adaptIntegrateSimplex(integrand_pf, S)
  # Integrand from PD = S_PPS(x2-x1)
  integrand_pd <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x))
    spps <- calc_surv(x, pps.ts$type, pps.ts$spec)
    vn*spps
    }
  int_pd <- stats::integrate(integrand_pd, 0, Tw)
  # Mean sojourn given starting vector
  soj <- starting[1] * int_pf$integral + starting[2] * int_pd$value
  return(soj)
}

#' Safely calculate restricted mean duration in progressed disease state for clock reset state transition model
#' @description Calculates the mean duration in the progressed disease state for the clock reset state transition model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over. Wrapper with 'possibly' of [rmd_pd_stm_cr]. This function is called by [calc_allrmds].
#' @param ... Pass-through to [rmd_pd_stm_cr]
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
prmd_pd_stm_cr <- purrr::possibly(rmd_pd_stm_cr, otherwise=NA_real_)

#' Restricted mean duration in progressed disease state for clock forward state transition model
#' @description Calculates the mean duration in the progressed disease state for the clock forward state transition model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over.
#' @inheritParams rmd_pf_stm
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
#' @seealso Used safely as [prmd_pd_stm_cf] by [calc_allrmds]
#' @export
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' # Find mean(s)
#' rmd_pd_stm_cf(dpam=params)
rmd_pd_stm_cf <- function(dpam, Ty=10, starting=c(1, 0, 0), lifetable=NA, discrate=0) {
  # Declare error if starting[2]>0 and lifetable!=NA
  if ((is.data.frame(lifetable))*starting[2]!=0) stop("Cannot apply life-table constraint when PD state not initially empty")
  # Declare local variables
  Tw <- ttp.ts <- ppd.ts <- pps.ts <- NULL
  S <- int_pf <- int_pd <- soj <- NULL
  # Bound to aid integration in weeks
  Tw <- convert_yrs2wks(Ty)
  # Normalize starting vector
  starting <- starting/sum(starting)
  # Pull out type and spec for TTP, PPD and PPS_CF
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ppd.ts <- convert_fit2spec(dpam$ppd)
  pps.ts <- convert_fit2spec(dpam$pps_cf)
  # Integrand PF = S_TTP(x1) * S_PPD(x1) * h_TTP(x1) * S_PPS_CF(x1, x2)
  # where S_PPS_CF = S_OS(x2)/S_OS(x1)
  # subject to lifetable maximum
  integrand_pf <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x[2]))
    sppd <- calc_surv(x[1], ppd.ts$type, ppd.ts$spec)
    fttp <- calc_dens(x[1], ttp.ts$type, ttp.ts$spec)
    sos1 <- calc_surv(x[1], pps.ts$type, pps.ts$spec)
    sos2 <- calc_surv(x[2], pps.ts$type, pps.ts$spec)
    spps <- sos2/sos1
    if (is.data.frame(lifetable)) {
      gens1 <- calc_ltsurv(convert_wks2yrs(x[1]), lifetable)
      gens2 <- calc_ltsurv(convert_wks2yrs(x[2]), lifetable)
      sppd <- pmin(sppd, gens1)
      spps <- pmin(spps, gens2/gens1)
      }
    if (sos1==0) 0 else {vn*sppd*fttp*spps}
  }
  S <- cbind(c(0,0),c(0, Tw),c(Tw, Tw))
  int_pf <- SimplicialCubature::adaptIntegrateSimplex(integrand_pf, S)
  # Integrand PD = S_OS(x2)/S_OS(x1)
  integrand_pd <- function(x) {
    calc_surv(x, pps.ts$type, pps.ts$spec)
  }
  int_pd <- stats::integrate(integrand_pd, 0, Tw)
  # Mean sojourn given starting vector
  soj <- starting[1] * int_pf$integral + starting[2] * int_pd$value
  return(soj)
}

#' Safely calculate restricted mean duration in progressed disease state for clock forward state transition model
#' @description Calculates the mean duration in the progressed disease state for the clock forward state transition model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over. Wrapper with 'possibly' of [rmd_pd_stm_cf]. This function is called by [calc_allrmds].
#' @param ... Pass-through to [rmd_pd_stm_cf]
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
prmd_pd_stm_cf <- purrr::possibly(rmd_pd_stm_cf, otherwise=NA_real_)

#' Restricted mean duration in progression free state for the partitioned survival model
#' @description Calculates the mean duration in the progression free state for the partitioned survival model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over.
#' @inheritParams rmd_pf_stm
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @seealso Used safely as [prmd_pf_psm] by [calc_allrmds]
#' @include basics.R
#' @export
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' # Find mean(s)
#' rmd_pf_psm(dpam=params)
rmd_pf_psm <- function(dpam, Ty=10, starting=c(1, 0, 0), lifetable=NA, discrate=0) {
  # Declare local variables
  Tw <- pfs.ts <- rmd <- NULL
  # Bound to aid integration in weeks
  Tw <- convert_yrs2wks(Ty)
  # Normalize starting vector
  starting <- starting/sum(starting)
  # Pull out type and spec for PFS
  pfs.ts <- convert_fit2spec(dpam$pfs)
  ttp.ts <- convert_fit2spec(dpam$ttp)
  # Create an integrand for PF survival
  integrand_pf <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x))
    pf_psm <- calc_surv(x, pfs.ts$type, pfs.ts$spec)
    if (is.data.frame(lifetable)) {
      ttp <- calc_surv(x, ttp.ts$type, ttp.ts$spec)
      osmax <- calc_ltsurv(convert_wks2yrs(x), lifetable)
      pf_psm <- pmin(pf_psm, ttp*osmax)
    }
    vn * pf_psm
  }
  int_pf <- stats::integrate(integrand_pf, 0, Tw)          
  # Finally multiply by starting population in PF
  starting[1] * int_pf$value
}

#' Safely calculate restricted mean duration in progression free state for the partitioned survival model
#' @description Calculates the mean duration in the progression free state for the partitioned survival model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over. Wrapper with 'possibly' of [rmd_pf_psm]. This function is called by [calc_allrmds].
#' @param ... Pass-through to [rmd_pf_psm]
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
prmd_pf_psm <- purrr::possibly(rmd_pf_psm, otherwise=NA_real_)

#' Restricted mean duration for overall survival in the partitioned survival model
#' @description Calculates the mean duration alive (i.e. in either progression free or progressed disease states) for the partitioned survival model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over.
#' @inheritParams rmd_pf_stm
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
#' @seealso Used safely as [prmd_os_psm] by [calc_allrmds]
#' @export
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' rmd_os_psm(params)
rmd_os_psm <- function(dpam, Ty=10, starting=c(1, 0, 0), lifetable=NA, discrate=0) {
  # Declare local variables
  Tw <- os.ts  <- NULL
  # Bound to aid integration in weeks
  Tw <- convert_yrs2wks(Ty)
  # Normalize starting vector
  starting <- starting/sum(starting)
  # Pull out type and spec for OS
  os.ts <- convert_fit2spec(dpam$os)
  # Create an integrand for overall survival
  integrand_os <- function(x) {
    vn <- (1+discrate)^(-convert_wks2yrs(x))
    os_psm <- calc_surv(x, os.ts$type, os.ts$spec)
    if (is.data.frame(lifetable)) {
      osmax <- calc_ltsurv(convert_wks2yrs(x), lifetable)
      os_psm <- pmin(os_psm, osmax)                     
    }
    vn*os_psm
  }
  int_os <- stats::integrate(integrand_os, 0, Tw)
  # Finally multiply by starting population in OS
  (starting[1] + starting[2]) * int_os$value
}

#' Safely calculate restricted mean duration for overall survival in the partitioned survival model
#' @description Calculates the mean duration alive (i.e. in either progression free or progressed disease states) for the partitioned survival model. Requires a carefully formatted list of fitted survival regressions for necessary endpoints, and the time duration to calculate over. Wrapper with 'possibly' of [rmd_os_psm]. This function is called by [calc_allrmds].
#' @param ... Pass-through to [rmd_os_psm]
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R
prmd_os_psm <- purrr::possibly(rmd_os_psm, otherwise=NA_real_)

#' Calculate restricted mean durations for each health state and all three models
#' @description Calculate restricted mean durations for each health state (progression free and progressed disease) for all three models (partitioned survival, clock forward state transition model, clock reset state transition model).
#' @param simdat Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
#' @param inclset Vector to indicate which patients to include in analysis
#' @param cuttime Time cutoff - this is nonzero for two-piece models.
#' @param Ty Time duration over which to calculate. Assumes input is in years, and patient-level data is recorded in weeks.
#' @param dpam List of statistical fits to each endpoint required in PSM, STM-CF and STM-CR models.
#' @param lifetable Optional, a life table. Columns must include `lttime` (time in years, or 52.18 times shorter than the time index elsewhere, starting from zero) and `lx`
#' @param discrate Discount rate (% per year)
#' @return List of detailed numeric results
#' - cutadj indicates the survival function and area under the curves for PFS and OS up to the cutpoint
#' - results provides results of the restricted means calculations, by model and state.
#' @include basics.R
#' @seealso Restricted means are provided by [rmd_pf_psm()], [rmd_os_psm()], [rmd_pf_stm()], [rmd_pd_stm_cf()] and [rmd_pd_stm_cr()]. The function [calc_allrmds_boot] provides a version for bootstrapping.
#' @export
#' @examples
#' # Create dataset and fit survival models (splines)
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' # Find mean(s)
#' calc_allrmds(bosonc, dpam=params)
calc_allrmds <- function(simdat,
                         inclset = 0,
                         dpam,
                         cuttime = 0,
                         Ty = 10,
                         lifetable = NA,
                         discrate = 0) {
  # Declare local variables
  ds <- ts.ppd <- fit.ppd <- ts.ttp <- fit.ttp <- NULL
  ts.pfs <- fit.pfs <- ts.os <- fit.os <- NULL
  ts.pps_cf <- fit.pps_cf <- ts.pps_cr <- fit.pps_cr <- NULL
  ds <- chvalid <- pf_km <- os_km <- NULL
  pfdat <- pfarea <- pfsurv <- osdat <- osarea <- ossurv <- NULL
  adjTy <- pf_psm_post <- os_psm_post <- pf_stm_post <- pd_stmcf_post <- pd_stmcr_post <- NULL
  pf_psm <- os_psm <- pf_stm <- os_stm_cf <- os_stm_cr <- NULL
  pd_psm <- pd_stm_cf <- pd_stm_cr <- NULL
  vec_pf <- vec_pd <- vec_os <- vec_model <- res <- cutadj <- NULL
  # Return vector for two-piece adjustment
  cutadj <- list(pf_area=pfarea, pf_surv=pfsurv, os_area=osarea, os_surv=ossurv)
  # Calculate rest through differences
  pd_psm_post <- os_psm_post-pf_psm_post
  os_stmcf_post <- pf_stm_post+pd_stmcf_post
  os_stmcr_post <- pf_stm_post+pd_stmcr_post
  # For a bootstrap sample, refit all distributions
  if (inclset[1]!=0) {
    simdat <- simdat[inclset,]
    ds <- create_extrafields(simdat, cuttime)
    dspps <- ds |> dplyr::filter(pps.durn>0, ttp.flag==1) 
    # Fit chosen distributions to each endpoint - PPD
    ts.ppd <- convert_fit2spec(dpam$ppd)
    fit.ppd <- fit_mods(durn1 = ds$tzero,
                       durn2 = ds$ppd.durn,
                       evflag = ds$ppd.flag,
                       type = ts.ppd$type,
                       spec = ts.ppd$spec
                       )[[1]]
    # TTP
    ts.ttp <- convert_fit2spec(dpam$ttp)
    fit.ttp <- fit_mods(durn1 = ds$tzero,
                       durn2 = ds$ttp.durn,
                       evflag = ds$ttp.flag,
                       type = ts.ttp$type,
                       spec = ts.ttp$spec
                       )[[1]]
    # PFS
    ts.pfs <- convert_fit2spec(dpam$pfs)
    fit.pfs <- fit_mods(durn1 = ds$tzero,
                       durn2 = ds$pfs.durn,
                       evflag = ds$pfs.flag,
                       type = ts.pfs$type,
                       spec = ts.pfs$spec
                       )[[1]]
    # OS
    ts.os <- convert_fit2spec(dpam$os)
    fit.os <- fit_mods(durn1 = ds$tzero,
                      durn2 = ds$os.durn,
                      evflag = ds$os.flag,
                      type = ts.os$type,
                      spec = ts.os$spec
                      )[[1]]
    # PPS CF - requires two time values
    ts.pps_cf <- convert_fit2spec(dpam$pps_cf)
    fit.pps_cf <- fit_mods(durn1 = dspps$ttp.durn,
                      durn2 = dspps$os.durn,
                      evflag = dspps$pps.flag,
                      type = ts.pps_cf$type,
                      spec = ts.pps_cf$spec
                      )[[1]]
    # PPS CR
    ts.pps_cr <- convert_fit2spec(dpam$pps_cr)
    fit.pps_cr <- fit_mods(durn1 = dspps$tzero,
                       durn2 = dspps$pps.durn,
                       evflag = dspps$pps.flag,
                       type = ts.pps_cr$type,
                       spec = ts.pps_cr$spec
                       )[[1]]
    # Bring together the best fits for each endpoint in a list
    dpam <- list(ppd=fit.ppd$result,
                   ttp=fit.ttp$result,
                   pfs=fit.pfs$result,
                   os=fit.os$result,
                   pps_cf=fit.pps_cf$result,
                   pps_cr=fit.pps_cr$result)
  } else {
    ds <- create_extrafields(simdat, cuttime)
  }
  # Check calculations valid
  chvalid <- is.na(dpam[1])==FALSE
  # Calculate mean duration and survival probability up to cutoff time
  if (cuttime>0) {
    pf_km <- survival::survfit(
      survival::Surv(pfs.odurn, pfs.flag) ~ 1, data=ds)
    os_km <- survival::survfit(
      survival::Surv(os.odurn, os.flag) ~ 1, data=ds)
    # PF calculations
    pfdat <- tidyr::as_tibble(cbind(time=pf_km$time, surv=pf_km$surv)) |>
      dplyr::mutate(
        row = dplyr::row_number(),
        itime = dplyr::if_else(.data$row==1, .data$time, .data$time-dplyr::lag(.data$time)),
        incl = dplyr::if_else(.data$time<cuttime, 1, 0),
        area = .data$incl*.data$surv*.data$itime
      )
    pfarea <- sum(pfdat$area)
    pfsurv <- min(pfdat[pfdat$incl==1,]$surv)
    # OS calculations
    osdat <- tidyr::as_tibble(cbind(time=os_km$time, surv=os_km$surv)) |>
      dplyr::mutate(
        row = dplyr::row_number(),
        itime = dplyr::if_else(.data$row==1, .data$time, .data$time-dplyr::lag(.data$time)),
        incl = dplyr::if_else(.data$time<cuttime, 1, 0),
        area = .data$incl*.data$surv*.data$itime
      )
    osarea <- sum(osdat$area)
    ossurv <- min(osdat[osdat$incl==1,]$surv)
  } else {
    # No adjustments if cuttime<=0
    pfarea <- osarea <-0
    pfsurv <- ossurv <-1
  }
  starting <- c(pfsurv, ossurv-pfsurv, 1-ossurv)
  if (chvalid) {
    # Call functions to calculate mean durations
    adjTy <- Ty - cuttime*7/365.25
    pf_psm_post <- prmd_pf_psm(dpam, Ty=adjTy, starting=starting, lifetable=lifetable, discrate=discrate)
    os_psm_post <- prmd_os_psm(dpam, Ty=adjTy, starting=starting, lifetable=lifetable, discrate=discrate)
    pf_stm_post <- prmd_pf_stm(dpam, Ty=adjTy, starting=starting, lifetable=lifetable, discrate=discrate)
    pd_stmcf_post <- prmd_pd_stm_cf(dpam, Ty=adjTy, starting=starting, lifetable=lifetable, discrate=discrate)
    pd_stmcr_post <- prmd_pd_stm_cr(dpam, Ty=adjTy, starting=starting, lifetable=lifetable, discrate=discrate)
  }
  # Calculating overall means
  pf_psm <- pfarea + pf_psm_post
  os_psm <- osarea + os_psm_post
  pf_stm <- pfarea + pf_stm_post
  os_stm_cf <- osarea + pf_stm_post + pd_stmcf_post
  os_stm_cr <- osarea + pf_stm_post + pd_stmcr_post
  pd_psm <- os_psm-pf_psm
  pd_stm_cf <- os_stm_cf-pf_stm
  pd_stm_cr <- os_stm_cr-pf_stm
  # Warning if PD LYs are zero
  if(pd_stm_cr<0.1 | is.na(pd_stm_cr)) {warning("Mean LY for PD state in STM_CR may not be calculable. Try a shorter time horizon.")}
  if(pd_stm_cf<0.1 | is.na(pd_stm_cf)) {warning("Mean LY for PD state in STM_CF may not be calculable. Try a shorter time horizon.")}
  # Return matrix of PF/PD/OS by method
  vec_pf <- c(pf_psm, pf_stm, pf_stm)
  vec_pd <- c(pd_psm, pd_stm_cf, pd_stm_cr)
  vec_os <- c(os_psm, os_stm_cf, os_stm_cr)
  vec_model <- c("PSM", "STM-CF", "STM-CR")
  res <- tidyr::tibble(model=vec_model, pf=vec_pf, pd=vec_pd, os=vec_os)
  # Return vector for two-piece adjustment
  cutadj <- list(pf_area=pfarea, pf_surv=pfsurv, os_area=osarea, os_surv=ossurv)
  return(list(cutadj=cutadj, results=res))
}

#' Wrapper to enable bootstrap sampling of restricted mean durations for each health state and all three models.
#' @description Wrapper function to [calc_allrmds] to enable bootstrap sampling of calculations of restricted mean durations for each health state (progression free and progressed disease) for all three models (partitioned survival, clock forward state transition model, clock reset state transition model).
#' @inheritParams calc_allrmds
#' @return Numeric vector of restricted mean durations - PF for each model (PSM, STM-CF, STM-CR), then PD, then OS.
#' @include basics.R
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#' )
#' calc_allrmds_boot(bosonc, dpam=params)
calc_allrmds_boot <- function(simdat,
                              inclset = 0,
                              dpam,
                              cuttime = 0,
                              Ty = 10,
                              lifetable = NA,
                              discrate = 0) {
  if (inclset[1]==0) {inclset <- 1:length(simdat$ptid)}
  mb <- calc_allrmds(simdat = simdat,
                     inclset = inclset,
                     dpam = dpam,
                     cuttime = cuttime,
                     Ty = Ty,
                     lifetable = lifetable,
                     discrate = discrate)
  c(mb$results$pf, mb$results$pd, mb$results$os)
}
