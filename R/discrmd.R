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
# Discretized restricted mean durations
# discrmd.R
# =====================================
#

#' Discretized Restricted Mean Duration calculation for Partitioned Survival Model
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a Partitioned Survival Model structure.
#' @param dpam List of survival regressions for model endpoints. These must include time to progression (TTP) and pre-progression death (PPD).
#' @param Ty Time duration over which to calculate (defaults to 10 years). Assumes input is in years, and patient-level data is recorded in weeks.
#' @param discrate Discount rate (%) per year (defaults to zero).
#' @param lifetable Optional. The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param timestep Optional, defaults to one (week).
#' @return List containing:
#' - pf: RMD in PF state
#' - pd: RMD in PD state
#' - os: RMD in either alive state
#' @seealso [drmd_stm_cf()] [drmd_stm_cr()]
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
#' drmd_psm(dpam=params)
#' # Add a lifetable constraint
#' ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
#' drmd_psm(dpam=params, lifetable=ltable)
drmd_psm <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- pfprob <- osprob <- adjosprob <- adjprob <- vn <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(1:floor(Tw/timestep)) + timestep/2
  # Membership probabilities with lifetable constraint
  pfprob <- prob_pf_psm(tvec, dpam)
  osprob <- prob_os_psm(tvec, dpam)
  adjosprob <- constrain_survprob(osprob, lifetable=lifetable, timevec=tvec)
  adjpfprob <- pfprob * adjosprob/osprob
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(adjpfprob*vn) * timestep
  os <- sum(adjosprob*vn) * timestep
  # Return values
  return(list(pf=pf, pd=os-pf, os=os))
}

#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Forward structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Forward structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_psm()] [drmd_stm_cr()]
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
#' drmd_stm_cf(dpam=params)
#' # Add a lifetable constraint
#' ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
#' drmd_stm_cf(dpam=params, lifetable=ltable)
drmd_stm_cf <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- ppd.ts <- ttp.ts <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(1:floor(Tw/timestep)) + timestep/2
  # Pull out type and spec for PPD and TTP
  ppd.ts <- convert_fit2spec(dpam$ppd)
  ttp.ts <- convert_fit2spec(dpam$ttp)
  # Calculate S_PPD, S_TTP and S_OS
  sppd <- tvec |> purrr::map_dbl(~calc_surv(.x, ppd.ts$type, ppd.ts$spec))
  sttp <- tvec |> purrr::map_dbl(~calc_surv(.x, ttp.ts$type, ttp.ts$spec))
  # Next line is the difference with STM-CR
  sos <- prob_os_stm_cf(tvec, dpam)
  # Apply constraints to S_PPD and S_OS
  adjsppd <- constrain_survprob(sppd, lifetable=lifetable, timevec=tvec)
  adjos <- constrain_survprob(sos, lifetable=lifetable, timevec=tvec)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(sttp*adjsppd*vn) * timestep
  os <- sum(adjos*vn) * timestep
  # Return values
  return(list(pf=pf, pd=os-pf, os=os))
}

#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Reset structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Reset structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_stm_cf()] [drmd_psm()]
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
#' drmd_stm_cr(dpam=params)
#' # Add a lifetable constraint
#' ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
#' drmd_stm_cr(dpam=params, lifetable=ltable)
drmd_stm_cr <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- ppd.ts <- ttp.ts <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(1:floor(Tw/timestep)) + timestep/2
  # Pull out type and spec for PPD and TTP
  ppd.ts <- convert_fit2spec(dpam$ppd)
  ttp.ts <- convert_fit2spec(dpam$ttp)
  # Calculate S_PPD, S_TTP and S_OS
  sppd <- tvec |> purrr::map_dbl(~calc_surv(.x, ppd.ts$type, ppd.ts$spec))
  sttp <- tvec |> purrr::map_dbl(~calc_surv(.x, ttp.ts$type, ttp.ts$spec))
  # Next line is the difference with STM-CF
  sos <- prob_os_stm_cr(tvec, dpam)
  # Apply constraints to S_PPD and S_OS
  adjsppd <- constrain_survprob(sppd, lifetable=lifetable, timevec=tvec)
  adjos <- constrain_survprob(sos, lifetable=lifetable, timevec=tvec)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(sttp*adjsppd*vn) * timestep
  os <- sum(adjos*vn) * timestep
  # Return values
  return(list(pf=pf, pd=os-pf, os=os))
}


