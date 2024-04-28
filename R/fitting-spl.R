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
# Functions involved in fitting splines
# fitting-spl.R
# ===========================================

#' Fit survival regressions with multiple spline regressions
#' @description Fits survival regressions with [flexsurv::flexsurvspline] to multiple Royston-Parmar splines specifications in one function call.
#' @param durn1 First time point, corresponds with time in [survival::Surv()]. For right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.
#' @param durn2 Second time point, corresponds with time2 in [survival::Surv()]. The ending time of the interval for interval censored or counting process data only.
#' @param evflag Event flag, corresponds with event in [survival::Surv]. The status indicator is normally 0=alive, 1=dead.
#' @param knot_set is a vector of the numbers of knots to consider, following [flexsurv::flexsurvspline()]).
#' @param scale_set is a vector of the spline scales to consider, following [flexsurv::flexsurvspline()]).
#' @param expvar Explanatory variable for modeling of PPS
#' @return A list by distribution, each containing two components:
#' - result: A list of class [flexsurv::flexsurvspline] containing information about the fitted model.
#' - error: Any error message returned on fitting the regression (NULL indicates no error).
#' @seealso [flexsurv::flexsurvspline] and [survival::Surv()]
#' @noRd
fit_mods_spl <- function(durn1, durn2=NA, evflag,
                         knot_set=1:3,
                         scale_set=c("hazard", "odds", "normal"),
                         expvar = NA
                         ) {
  # Declare local variables
  multispl <- entry <- NULL
  # Return nothing if no knots or scales specified
  if (sum(is.na(knot_set))!=0) {return(NA)}
  if (sum(is.na(scale_set))!=0) {return(NA)}
  # Safely run spline fittings (i.e. capture when regression fails)
  sflexsurvspline <- purrr::safely(.f = flexsurv::flexsurvspline)
  # Create an empty list
  multispl <- vector("list", length(knot_set) * length(scale_set))
  # Fill list entries with spline models
  for (k in seq(knot_set)) {
    for (s in seq(scale_set)) {
      entry <- s+(k-1)*length(knot_set)
      # If there is an explanatory variable
      if (!is.na(expvar[1])) {
        multispl[[entry]] <- sflexsurvspline(
            formula = survival::Surv(time = durn1,
                                 time2 = durn2,
                                 event = evflag) ~ {{ expvar }},
            k = knot_set[k],
            scale = scale_set[s]
        ) } else {
          # Or if there is not
          multispl[[entry]] <- sflexsurvspline(
            formula = survival::Surv(time = durn1,
                                     time2 = durn2,
                                     event = evflag) ~ 1,
            k = knot_set[k],
            scale = scale_set[s]
          ) }
    }
  }
  return(multispl)
}  

#' Fit multiple spline regressions to the multiple required endpoints
#' @description Fits multiple survival regressions, according to the distributions stipulated, to the multiple endpoints required in fitting partitioned survival analysis, clock forward and clock reset semi-markov models.
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
#' @param knot_set is a vector of the numbers of knots to consider, following [flexsurv::flexsurvspline()]).
#' @param scale_set is a vector of the spline scales to consider, following [flexsurv::flexsurvspline()]).
#' @param expvar Explanatory variable for modeling of PPS
#' @return A list by endpoint, then distribution, each containing two components:
#' - result: A list of class [flexsurv::flexsurvspline] containing information about the fitted model.
#' - error: Any error message returned on fitting the regression (NULL indicates no error).
#' Also, the given cuttime.
#' @export
#' @seealso Parametric modeling is handled by [fit_ends_mods_par()]
#' @examples
#' \donttest{
#' # Create dataset in suitable form using bos dataset from the flexsurv package
#' bosonc <- create_dummydata("flexbosms")
#' fit_ends_mods_spl(bosonc, expvar=bosonc$ttp.durn)
#' }
fit_ends_mods_spl <- function(simdat,
                              knot_set=1:3,
                              scale_set=c("hazard", "odds", "normal"),
                              expvar = NA
                              ) {
  # Declare local variables
  ds <- dspps <- pps.durn <- fits.ppd <- fits.ttp <- fits.pfs <- NULL
  fits.os <- fits.pps_cf <- fits.pps_cr <- NULL
  # Derive additional fields, as with regular function
  ds <- create_extrafields(simdat, cuttime=0)
  # For PPS analysis, require there to be a known progression event, plus a positive PPS
  dspps <- ds |> dplyr::filter(.data$pps.durn>0, .data$ttp.flag==1)
  # Captures lists of fitted models to each endpoint
  fits.ppd <- fit_mods_spl(durn1 = ds$tzero,
                       durn2 = ds$ppd.durn,
                       evflag = ds$ppd.flag,
                       knot_set = knot_set,
                       scale_set = scale_set
                       )
  fits.ttp <- fit_mods_spl(durn1 = ds$tzero,
                       durn2 = ds$ttp.durn,
                       evflag = ds$ttp.flag,
                       knot_set = knot_set,
                       scale_set = scale_set
                        )
  fits.pfs <- fit_mods_spl(durn1 = ds$tzero,
                       durn2 = ds$pfs.durn,
                       evflag = ds$pfs.flag,
                       knot_set = knot_set,
                       scale_set = scale_set
                        )
  fits.os <- fit_mods_spl(durn1 = ds$tzero,
                      durn2 = ds$os.durn,
                      evflag = ds$os.flag,
                      knot_set = knot_set,
                      scale_set = scale_set
                        )
  # CR requires two time values
  fits.pps_cf <- fit_mods_spl(durn1=dspps$ttp.durn,
                          durn2=dspps$os.durn,
                          evflag = dspps$pps.flag,
                          knot_set = knot_set,
                          scale_set = scale_set,
                          expvar = expvar
                          )
  fits.pps_cr <- fit_mods_spl(durn1 = dspps$tzero,
                      durn2 = dspps$pps.durn,
                      evflag = dspps$pps.flag,
                      knot_set = knot_set,
                      scale_set = scale_set,
                      expvar = expvar
                        )
  # Return a list of all the fits
  list(ttp=fits.ttp,
       ppd=fits.ppd,
       pfs=fits.pfs,
       os=fits.os,
       pps_cf=fits.pps_cf,
       pps_cr=fits.pps_cr)
}
