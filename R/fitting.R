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
# Functions involved in fitting parametric distributions
# fitting.R
# ==================================================================
#
#' Checks whether the Hessian for a survival regression is positive-definite
#' @description Checks whether the Hessian matrix returned in a list after fitting a survival regression with flexsurvreg is positive-definite.
#' @param fitlist is a list returned after running [flexsurv::flexsurvreg()] describing a fitted survival model
#' @return logical: TRUE if Hessian matrix is positive definite, FALSE if not.
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_par(bosonc)
#' check_posdef(fits$pfs[[2]]$result)
check_posdef <- function(fitlist) {
  out <- tryCatch({
    det(chol(fitlist$opt$hessian))>0
    },
  error = function(cond) {
    message(paste("Hessian is not positive definite when fitted to the",
                  fitlist$result$dlist$name,
                  "distribution. \n"
                  ))
    FALSE
  },
  warning = function(cond) {FALSE},
  finally = {}
  )
  return(out)
}

#' Fit survival regressions with multiple parametric distributions
#' @description Fits survival regressions with [flexsurv::flexsurvreg] to multiple statistical distributions in one function call.
#' @param durn1 First time point, corresponds with time in [survival::Surv()]. For right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.
#' @param durn2 Second time point, corresponds with time2 in [survival::Surv()]. The ending time of the interval for interval censored or counting process data only.
#' @param evflag Event flag, corresponds with event in [survival::Surv]. The status indicator is normally 0=alive, 1=dead.
#' @param dists is a vector of all the distributions (named according to [flexsurv::flexsurvreg()]) for which statistical fits are requested.
#' @param expvar Explanatory variable for modeling of PPS
#' @return A list by distribution, each containing two components:
#' - result: A list of class [flexsurv::flexsurvreg] containing information about the fitted model.
#' - error: Any error message returned on fitting the regression (NULL indicates no error).
fit_mods_par <- function(durn1, durn2=NA, evflag, dists, expvar=NA) {
  # Return nothing if no distributions specified
  if (sum(is.na(dists))!=0) {return(NA)}
  # Safely run regressions (i.e. capture when regression fails)
  sflexsurvreg <- purrr::safely(.f = flexsurv::flexsurvreg)
  # Run regression for each distribution in dists
  # ... If there is an explanatory variable
  if (!is.na(expvar[1])) {
    dists |> purrr::map(~sflexsurvreg(
        formula = survival::Surv(time = durn1,
                             time2 = durn2,
                             event = evflag) ~ {{ expvar }},
        dist=.x)
        )
  } else {
    # ... Or if there is not
    dists |> purrr::map(~sflexsurvreg(
      formula = survival::Surv(time = durn1,
                               time2 = durn2,
                               event = evflag) ~ 1,
      dist=.x)
    )
  }
}

#' Fit survival regressions with multiple splines or parametric distributions
#' @param durn1 Start time
#' @param durn2 End time
#' @param evflag Event flag
#' @param type Type of model ("spl" for spline or "par" for parametric)
#' @param spec Specification of model
#' @return A list by distribution, each containing two components:
#' - result: A list containing information about the fitted model.
#' - error: Any error message returned on fitting the regression (NULL indicates no error).
fit_mods <- function(durn1, durn2=NA, evflag, type, spec) {
  if (type=="spl"){
    fit_mods_spl(durn1, durn2, evflag, knot_set=spec$k, scale_set=spec$scale)
  } else if (type=="par") {
    fit_mods_par(durn1, durn2, evflag, dists=spec$dist)
  } else {stop("Incorrect model type specified. Must be par (parameteric) or spl (splines).")}
}

#' Create the additional time-to-event endpoints, adjusting for cutpoint
#' @param ds Patient-level dataset
#' @param cuttime Time cutpoint
#' @return Dataset of complete patient-level dataset, adjusted for cutpoint
#' @export
#' @importFrom rlang .data
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' create_extrafields(bosonc, cuttime=10)
create_extrafields <- function(ds, cuttime=0){
  ds <- ds |>
    dplyr::rename(
      ttp.odurn = "ttp.durn",
      pfs.odurn = "pfs.durn",
      os.odurn = "os.durn"
    ) |>
    dplyr::mutate(
      # Calculate time after time cut-off
      ttp.durn = pmax(0, .data$ttp.odurn - cuttime),
      pfs.durn = pmax(0, .data$pfs.odurn - cuttime),
      os.durn = pmax(0, .data$os.odurn - cuttime),
      tzero = 0,
      # Calculate other endpoints
      pps.odurn = pmax(.data$os.odurn - .data$ttp.odurn, 0),
      pps.durn = pmax(.data$os.durn - .data$ttp.durn, 0),
      pps.flag = .data$ttp.flag * .data$os.flag,
      ppd.durn = .data$ttp.durn,
      ppd.flag = (1 - .data$ttp.flag) * .data$pfs.flag
      )
  return(ds)
}

#' Fit multiple parametric survival regressions to the multiple required endpoints
#' @description Fits multiple parametric survival regressions, according to the distributions stipulated, to the multiple endpoints required in fitting partitioned survival analysis, clock forward and clock reset semi-markov models.
#' @param simdat Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
#'
#' Survival data for all other endpoints (time to progression, pre-progression death, post-progression survival) are derived from PFS and OS.
#' @param cuttime Cut-off time for a two-piece model, equals zero for one-piece models.
#' @param ppd.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Pre-Progression Death (PPD).
#' @param ttp.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Time To Progression (TTP).
#' @param pfs.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Progression-Free Survival (PFS).
#' @param os.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Overall Survival (OS).
#' @param pps_cf.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Post Progression Survival, where time is from baseline (clock forward).
#' @param pps_cr.dist Vector of distributions (named per [flexsurv::flexsurvreg()]) to be fitted to Post Progression Survival, where time is from progression (clock reset).
#' @param expvar Explanatory variable for modeling of PPS
#' @return A list by endpoint, then distribution, each containing two components:
#' - result: A list of class *flexsurvreg* containing information about the fitted model.
#' - error: Any error message returned on fitting the regression (NULL indicates no error).
#' @export
#' @seealso Spline modeling is handled by [fit_ends_mods_spl()]
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fit_ends_mods_par(bosonc, expvar=bosonc$ttp.durn)
fit_ends_mods_par <- function(simdat,
                cuttime = 0,
                ppd.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                ttp.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                pfs.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                os.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                pps_cf.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                pps_cr.dist = c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz"),
                expvar = NA) {
  # Declare local variables
  ds <- dspps <- pps.durn <- NULL
  fits.ppd <- fits.ttp <- fits.pfs <- fits.os <- fits.pps_cf <- fits.pps_cr <- NULL
  # Derive additional fields, as with regular function
  ds <- create_extrafields(simdat, cuttime)
  # For PPS analysis, require there to be a known progression event, plus a positive PPS
  dspps <- ds |> dplyr::filter(.data$pps.durn>0, .data$ttp.flag==1)
  # Captures lists of fitted models to each endpoint
  fits.ppd <- fit_mods_par(durn1 = ds$tzero,
                       durn2 = ds$ppd.durn,
                       evflag = ds$ppd.flag,
                       dists = ppd.dist)
  fits.ttp <- fit_mods_par(durn1 = ds$tzero,
                       durn2 = ds$ttp.durn,
                       evflag = ds$ttp.flag,
                       dists = ttp.dist)
  fits.pfs <- fit_mods_par(durn1 = ds$tzero,
                       durn2 = ds$pfs.durn,
                       evflag = ds$pfs.flag,
                       dists = pfs.dist)
  fits.os <- fit_mods_par(durn1 = ds$tzero,
                      durn2 = ds$os.durn,
                      evflag = ds$os.flag,
                      dists = os.dist)
  # CR requires two time values
  fits.pps_cf <- fit_mods_par(durn1=dspps$ttp.durn,
                          durn2=dspps$os.durn,
                          evflag = dspps$pps.flag,
                          dists = pps_cf.dist,
                          expvar = expvar)
  fits.pps_cr <- fit_mods_par(durn1 = dspps$tzero,
                      durn2 = dspps$pps.durn,
                      evflag = dspps$pps.flag,
                      dists = pps_cr.dist,
                      expvar = expvar)
  # Return a list of all the fits
  list(ttp=fits.ttp,
       ppd=fits.ppd,
       pfs=fits.pfs,
       os=fits.os,
       pps_cf=fits.pps_cf,
       pps_cr=fits.pps_cr
       )
}

#' Find the "best" survival regression from a list
#' @description When there are multiple survival regressions fitted to the same endpoint and dataset, it is necessary to identify the preferred model. This function reviews the fitted regressions and selects that with the minimum Akaike or Bayesian Information Criterion (AIC, BIC), depending on user choice.
#' @param reglist List of fitted survival regressions to an endpoint and dataset.
#' @param crit Criterion to be used in selection of best fit, either "aic" (Akaike Information Criterion) or "bic" (Bayesian Information Criterion).
#' @return List of the single survival regression with the best fit.
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' # Pick some distributions to fit to all endpoints
#' dists <- c("exp", "llogis", "lnorm")
#' # Fit all distributions to all endpoints
#' fits <- fit_ends_mods_par(bosonc,
#'     cuttime = 0,
#'     ppd.dist = dists,
#'     ttp.dist = dists,
#'     pfs.dist = dists,
#'     os.dist = dists,
#'     pps_cr.dist = dists,
#'     pps_cf.dist = dists)
#' find_bestfit_par(fits$ttp, "aic")
#' find_bestfit_par(fits$os, "bic")
find_bestfit_par <- function(reglist, crit) {
  # Declare local variables
  noreg <- valid <- remain <- NULL
  npts <- pars <- aic <- loglik <- bic <- NULL
  ic <- chosen <- posdef <- conv <- dists <- NULL
  restab <- othtab <- NULL
  # Pick out the valid regressions (where valid==TRUE)
  noreg <- length(reglist)
  valid <- seq(noreg) |> purrr::map_lgl(~is.null(reglist[[.x]]$error))
  remain <- seq(noreg)*valid
  remain <- remain[remain>0]
  # Pull useful data from fitted regressions
  npts <- remain |>
    purrr::map_dbl(~reglist[[.x]]$result$N)
  pars <- remain |>
    purrr::map_dbl(~reglist[[.x]]$result$npars)
  aic <- remain |>
    purrr::map_dbl(~reglist[[.x]]$result$AIC)
  loglik <- remain |>
    purrr::map_dbl(~reglist[[.x]]$result$loglik)
  bic <- pars*log(npts)-2*loglik
  # Find min AIC or BIC depending on crit
  ic <- if (crit=="BIC"|crit=="bic") bic
        else if (crit=="AIC"|crit=="aic") aic
        else rep(NA, noreg)
  chosen <- which.min(ic)
  # Create a useful table
  posdef <- remain |>
    purrr::map_lgl(~check_posdef(reglist[[.x]]$result))
  conv <- remain |>
    purrr::map_lgl(~reglist[[.x]]$result$opt$convergence==0)
  dists <- remain |>
    purrr::map_chr(~reglist[[.x]]$result$dlist$name)
  restab <- tibble::tibble(id=remain, dists, npts, pars, loglik, conv, posdef, aic, bic)
  othtab <- tibble::tibble(id=seq(noreg)[valid==FALSE], dists="unknown", npts=0, pars=0, loglik=-Inf, conv=FALSE, posdef=FALSE, aic=Inf, bic=Inf)
  restab <- dplyr::add_row(restab, othtab) |>
    dplyr::mutate(
      rankaic = rank(aic),
      rankbic = rank(bic)
    )
  # Pull out just the result
  list(fit=reglist[[remain[chosen]]]$result, results=restab)
}
