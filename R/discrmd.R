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
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
#'
#' Survival data for all other endpoints (time to progression, pre-progression death, post-progression survival) are derived from PFS and OS.
#' @param dpam List of survival regressions for model endpoints. These must include time to progression (TTP) and pre-progression death (PPD).
#' @param psmtype Either "simple" or "complex" PSM formulation
#' @param Ty Time duration over which to calculate (defaults to 10 years). Assumes input is in years, and patient-level data is recorded in weeks.
#' @param discrate Discount rate (%) per year (defaults to zero).
#' @param lifetable Optional. The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param timestep Optional, defaults to one (week).
#' @return List containing:
#' - pf: RMD in PF state
#' - pd: RMD in PD state
#' - os: RMD in either alive state
#' @seealso [drmd_stm_cf()] [drmd_stm_cr()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_psm(ptdata=bosonc, dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_psm(ptdata=bosonc, dpam=params, lifetable=ltable)
drmd_psm <- function(ptdata, dpam, psmtype="simple", Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- pfprob <- osprob <- adjosprob <- adjfac <- adjprob <- vn <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- ceiling(convert_yrs2wks(Ty)/timestep)
  # Create dataset, starting with time vector, with half-cycle addition, and unconstrained probs
  ds <- tibble::tibble(
    tzero = timestep*(0:Tw),
    tmid = tzero + timestep/2,
    pfprob = prob_pf_psm(tzero, dpam),
    osprob = prob_os_psm(tzero, dpam),
    u_pf = pfprob,
    u_pd = osprob-pfprob
  )
  # Obtain all the hazards
  allh <- calc_haz_psm(timevar=ds$tmid, ptdata=ptdata, dpam=dpam, psmtype=psmtype)$adj
  # Derive the unconstrained PPD mortality probability
  ds$q_ppd <- 1-exp(-allh$ppd)
  # Derive the constrained life table
  ds$clx <- calc_ltsurv(convert_wks2yrs(ds$tzero), lifetable)
  # Other calculations on the dataset
  ds <- ds |>
    dplyr::mutate(
      # Derive the background mortality for this timepoint
      cqx = 1 - dplyr::lead(clx)/clx,
      # Derive the TTP probability (balancing item)
      q_pfs = 1 - dplyr::lead(u_pf)/u_pf,
      q_ttp = q_pfs-q_ppd,
      d_pf = u_pf * q_ppd,
      c_qpfs = q_ttp + pmax(q_ppd, cqx),
      # Derive the PPS mortality probability
      d_pfpd = u_pf + u_pd - dplyr::lead(u_pf) - dplyr::lead(u_pd),
      d_pps = d_pfpd - d_pf,
      q_pps = dplyr::if_else(u_pd==0, 0, d_pps / u_pd),
      # Constrained probabilities
      cqpfs = q_ttp + pmax(q_ppd, cqx),
      cqpps = pmax(q_pps, cqx),
      # Derive the constrained PF and PD memberships
      c_pf = u_pf,
      c_pd = u_pd,
    )
  # Derive the constrained PF and PD memberships
  for (t in 2:(Tw)) {
      ds$c_pf[t] = ds$c_pf[t-1] * (1-ds$cqpfs[t-1])
      ds$c_pd[t] = ds$c_pf[t-1] * ds$q_ttp[t-1] + ds$c_pd[t-1] * (1 - ds$cqpps[t-1])
  }
  # The final membership probabilities are zero
  ds$c_pf[Tw+1] <- ds$c_pd[Tw+1] <- 0
  # Final calculations on the dataset
  ds <- ds |>
      dplyr::mutate(
        # Discount factor
        vn = (1+discrate)^(-convert_wks2yrs(tmid)),
        # RMD components in each timestep
        rmd_pf = (c_pf + dplyr::lead(c_pf))/2 * vn * timestep,
        rmd_pd = (c_pd + dplyr::lead(c_pd))/2 * vn * timestep
      )
  ds$rmd_pf[Tw+1] <- ds$rmd_pd[Tw+1] <- 0
  # Calculate RMDs
  pf <- sum(ds$rmd_pf)
  pd <- sum(ds$rmd_pd)
  # Return values
  return(list(pf=pf, pd=pd, os=pf+pd, calc=ds))
}


#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Forward structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Forward structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_psm()] [drmd_stm_cr()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_stm_cf(dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_stm_cf(dpam=params, lifetable=ltable)
drmd_stm_cf <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- ceiling(convert_yrs2wks(Ty)/timestep)
  # Create dataset, starting with time vector, with half-cycle addition, and unconstrained probs
  ds <- tibble::tibble(
    tzero = timestep*(0:Tw),
    tmid = tzero + timestep/2,
    u_pf = prob_pf_stm(tzero, dpam),
    u_pd = prob_pd_stm_cf(tzero, dpam),
    # Calculate PPD hazard and probability
    h_ppd = calc_haz(tmid, survobj=dpam$ppd),
    q_ppd = 1-exp(-h_ppd),
    # Derive the constrained life table
    clx = calc_ltsurv(convert_wks2yrs(tzero), lifetable),
    # Derive the background mortality for this timepoint
    cqx = 1 - dplyr::lead(clx)/clx,
    # Derive the TTP probability (balancing item for PFS)
    q_pfs = 1 - dplyr::lead(u_pf)/u_pf,
    q_ttp = q_pfs-q_ppd,
    d_pf = u_pf * q_ppd,
    # Derive the PPS mortality probability
    d_pfpd = u_pf + u_pd - dplyr::lead(u_pf) - dplyr::lead(u_pd),
    d_pps = d_pfpd - d_pf,
    q_pps = dplyr::if_else(u_pd==0, 0, d_pps / u_pd),
    # Constrained probabilities
    cqpfs = q_ttp + pmax(q_ppd, cqx),
    cqpps = pmax(q_pps, cqx),
    # Initial constrained PF and PD
    c_pf = u_pf,
    c_pd = u_pd
  )
  # Derive the constrained PF and PD memberships
  for (t in 2:(Tw)) {
    ds$c_pf[t] = ds$c_pf[t-1] * (1-ds$cqpfs[t-1])
    ds$c_pd[t] = ds$c_pf[t-1] * ds$q_ttp[t-1] + ds$c_pd[t-1] * (1 - ds$cqpps[t-1])
  }
  # The final membership probabilities are zero
  ds$c_pf[Tw+1] <- ds$c_pd[Tw+1] <- 0
  # Final calculations on the dataset
  ds <- ds |>
    dplyr::mutate(
      # Discount factor
      vn = (1+discrate)^(-convert_wks2yrs(tmid)),
      # RMD components in each timestep
      rmd_pf = (c_pf + dplyr::lead(c_pf))/2 * vn * timestep,
      rmd_pd = (c_pd + dplyr::lead(c_pd))/2 * vn * timestep
    )
  ds$rmd_pf[Tw+1] <- ds$rmd_pd[Tw+1] <- 0
  # Calculate RMDs
  pf <- sum(ds$rmd_pf)
  pd <- sum(ds$rmd_pd)
  # Return values
  return(list(pf=pf, pd=pd, os=pf+pd, calc=ds))
}

#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Reset structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Reset structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_stm_cf()] [drmd_psm()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_stm_cr(dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_stm_cr(dpam=params, lifetable=ltable)
drmd_stm_cr <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- ceiling(convert_yrs2wks(Ty)/timestep)
  # Create dataset, starting with time vector, with half-cycle addition, and unconstrained probs
  ds <- tibble::tibble(
    tzero = timestep*(0:Tw),
    tmid = tzero + timestep/2,
    u_pf = prob_pf_stm(tzero, dpam),
    u_pd = prob_pd_stm_cf(tzero, dpam)
  )
  # Keep calculating membership probabilities
  ds <- ds |>
    dplyr::mutate(
    # Calculate PPD hazard and probability
    h_ppd = calc_haz(tmid, survobj=dpam$ppd),
    q_ppd = 1-exp(-h_ppd),
    # Derive the constrained life table
    clx = calc_ltsurv(convert_wks2yrs(tzero), lifetable),
    cqx = 1 - dplyr::lead(clx)/clx,
    # Derive the TTP probability (balancing item for PF)
    q_pfs = 1 - dplyr::lead(u_pf)/u_pf,
    q_ttp = q_pfs-q_ppd,
    d_pf = u_pf * q_ppd,
    # Derive the PPS mortality probability
    d_pfpd = u_pf + u_pd - dplyr::lead(u_pf) - dplyr::lead(u_pd),
    d_pps = d_pfpd - d_pf,
    q_pps = dplyr::if_else(u_pd==0, 0, d_pps / u_pd),
    # Constrained probabilities
    cqpfs = q_ttp + pmax(q_ppd, cqx),
    cqpps = pmax(q_pps, cqx),
    # Initial constrained PF and PD
    c_pf = u_pf,
    c_pd = u_pd
  )
  # Derive the constrained PF and PD memberships
  for (t in 2:(Tw)) {
    ds$c_pf[t] = ds$c_pf[t-1] * (1-ds$cqpfs[t-1])
    ds$c_pd[t] = ds$c_pf[t-1] * ds$q_ttp[t-1] + ds$c_pd[t-1] * (1 - ds$cqpps[t-1])
  }
  # The final membership probabilities are zero
  ds$c_pf[Tw+1] <- ds$c_pd[Tw+1] <- 0
  # Final calculations on the dataset
  ds <- ds |>
    dplyr::mutate(
      # Discount factor
      vn = (1+discrate)^(-convert_wks2yrs(tmid)),
      # RMD components in each timestep
      rmd_pf = (c_pf + dplyr::lead(c_pf))/2 * vn * timestep,
      rmd_pd = (c_pd + dplyr::lead(c_pd))/2 * vn * timestep
    )
  ds$rmd_pf[Tw+1] <- ds$rmd_pd[Tw+1] <- 0
  # Calculate RMDs
  pf <- sum(ds$rmd_pf)
  pd <- sum(ds$rmd_pd)
  # Return values
  return(list(pf=pf, pd=pd, os=pf+pd, calc=ds))
}
