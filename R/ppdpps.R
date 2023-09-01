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
# Functions relating to calculating pre and post progression death
# ppdpps.R
# ==================================================================

#' Derive pre and post-progression hazards of death under PSM
#' @description Derive the hazards of death pre- and post-progression under either simple or complex PSM formulations.
#' @param timevar Vector of times at which to calculate the hazards
#' @param ptdata Patient-level dataset
#' @param dpam List of fitted survival models for each endpoint
#' @param type Either "simple" or "complex" PSM formulation
#' @return List of pre, the pre-progression hazard, and post, the post-progression hazard
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#'   )
#' calc_haz_psm(0:10, ptdata=bosonc, dpam=params, type="simple")
#' calc_haz_psm(0:10, ptdata=bosonc, dpam=params, type="complex")
calc_haz_psm <- function(timevar, ptdata, dpam, type) {
  # PFS
  pfs.ts <- convert_fit2spec(dpam$pfs)
  pfs.type <- pfs.ts$type
  pfs.spec <- pfs.ts$spec
  # OS
  os.ts <- convert_fit2spec(dpam$os)
  os.type <- os.ts$type
  os.spec <- os.ts$spec
  # TTP
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ttp.type <- ttp.ts$type
  ttp.spec <- ttp.ts$spec
  # Events in dataset
  ne_pfs <- sum(ptdata$pfs.flag)
  ne_ttp <- sum(ptdata$ttp.flag)
  progfrac <- max(0, min(1, ne_ttp/ne_pfs))
  # Hazards and survival functions
  http <- calc_haz(timevar, ttp.type, ttp.spec)
  hpf <- calc_haz(timevar, pfs.type, pfs.spec)
  hos <- calc_haz(timevar, os.type, os.spec)
  sos <- calc_surv(timevar, os.type, os.spec)
  spf <- calc_surv(timevar, pfs.type, pfs.spec)
  # hppd varies by type
  hppd_simple <- pmax(0, pmin((1-progfrac)*hpf, sos*hos/spf))
  hppd_complex <- pmax(0, pmin(hpf-http, sos*hos/spf))
  hppd <- hppd_simple*(type=="simple") + hppd_complex*(type!="simple")
  hpps <- pmax(0, pmin((sos*hos-spf*hppd)/(sos-spf), 5000))
  # hpps and others
  hdiff <- hos-(spf*hppd+(sos-spf)*hpps)
  hpps[timevar==0] <- 0
  return(list(pre=hppd, post=hpps, diff=hdiff))
}

#' Derive PPS survival function under a PSM
#' @description Derive the PPS survival function under the simple or complex PSM formulation.
#' @param totime Vector of times to which the survival function is calculated
#' @param fromtime Vector of times from which the survival function is calculated
#' @param ptdata Patient-level dataset
#' @param dpam List of fitted survival models for each endpoint
#' @param type Either "simple" or "complex" PSM formulation
#' @return Vector of PPS survival function values
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#'   os = find_bestfit_spl(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
#'   )
#' calc_surv_psmpps(totime=1:10,
#'     fromtime=rep(1,10),
#'     ptdata=bosonc,
#'     dpam=params,
#'     type="simple")
calc_surv_psmpps <- function(totime, fromtime=0, ptdata, dpam, type="simple") {
  # Hazard function
  hazfn <- function(x) {
    calc_haz_psm(timevar=x, ptdata=ptdata, dpam=dpam, type=type)$post
  }
  # Cum hazard function
  cumhazfn <- function(ptrow) {
    if (fromtime[ptrow]>=totime[ptrow]) {return(0)}
    intres <- stats::integrate(hazfn,
              lower=fromtime[ptrow],
              upper=totime[ptrow]
              )
    if (intres$message=="OK") intres$value else NA
  }
  # Cumulative hazard for all patients
  cumH <- seq(fromtime) |> purrr::map_dbl(~cumhazfn(.x))
  # Survival
  exp(-cumH)
}
