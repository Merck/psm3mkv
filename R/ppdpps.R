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
#' @param dpam List of survival regressions for each endpoint:
#' - pre-progression death (PPD)
#' - time to progression (TTP)
#' - progression-free survival (PFS)
#' - overall survival (OS)
#' - post-progression survival clock forward (PPS-CF) and
#' - post-progression survival clock reset (PPS-CR).
#' @param psmtype Either "simple" or "complex" PSM formulation
#' @return List of pre, the pre-progression hazard, and post, the post-progression hazard
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
#'   )
#' calc_haz_psm(0:10, ptdata=bosonc, dpam=params, psmtype="simple")
#' calc_haz_psm(0:10, ptdata=bosonc, dpam=params, psmtype="complex")
calc_haz_psm <- function(timevar, ptdata, dpam, psmtype) {
  # Declare local variables
  pfs.ts <- pfs.type <- pfs.spec <- NULL
  os.ts <- os.type <- os.spec <- NULL
  ttp.ts <- ttp.type <- ttp.spec <- NULL
  ne_pfs <- ne_ttp <- progfrac <- NULL
  http <- hpf <- hos <- sos <- spf <- NULL
  hppd_simple <- hppd_complex <- hppd <- hpps <- hdiff <- NULL
  # PFS
  pfs.ts <- convert_fit2spec(dpam$pfs)
  pfs.type <- pfs.ts$type
  pfs.spec <- pfs.ts$spec
  hpf <- calc_haz(timevar, pfs.type, pfs.spec)
  spf <- calc_surv(timevar, pfs.type, pfs.spec)
  # OS
  os.ts <- convert_fit2spec(dpam$os)
  os.type <- os.ts$type
  os.spec <- os.ts$spec
  hos <- calc_haz(timevar, os.type, os.spec)
  sos <- calc_surv(timevar, os.type, os.spec)
  # TTP complex
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ttp.type <- ttp.ts$type
  ttp.spec <- ttp.ts$spec
  http_complex <- calc_haz(timevar, ttp.type, ttp.spec)
  # TTP simple
  ne_pfs <- sum(ptdata$pfs.flag)
  ne_ttp <- sum(ptdata$ttp.flag)
  progfrac <- max(0, min(1, ne_ttp/ne_pfs))
  http_simple <- progfrac*hpf
  # TTP
  typeflag <- ifelse(psmtype=="simple", 1, 0)
  http <- http_simple*typeflag + http_complex*(1-typeflag)
  # PPD
  hppd_unadj <- hpf-http
  hppd_simple <- pmax(0, pmin((1-progfrac)*hpf, sos*hos/spf))
  hppd_complex <- pmax(0, pmin(hpf-http, sos*hos/spf))
  hppd <- hppd_simple*typeflag + hppd_complex*(1-typeflag)
  # PPS
  hpps_unadj <- (sos*hos-spf*hppd)/(sos-spf)
  hpps <- pmax(0, pmin(hpps_unadj, 5000))
  hpps[timevar==0] <- 0
  # Diff
  # hdiff <- hos-(spf*hppd+(sos-spf)*hpps)
  # Adjusted for caps and collars
  hadj <- list(
    ttp = http,
    ppd = hppd,
    pfs = http+hppd,
    os = spf*hppd + (sos-spf)*hpps,
    pps = hpps
  )
  # Unadjusted for caps and collars
  hunadj <- list(
    ttp = http,
    ppd = hpf-http,
    pfs = hpf,
    os = hos,
    pps = hpps_unadj
  )
  return(list(adj=hadj, unadj=hunadj))
}

#' Derive PPS survival function under a PSM
#' @description Derive the post-progression survival (PPS) function under the simple or complex PSM formulation.
#' @param totime Vector of times to which the survival function is calculated
#' @param fromtime Vector of times from which the survival function is calculated
#' @param ptdata Patient-level dataset
#' @param dpam List of fitted survival models for each endpoint
#' @param psmtype Either "simple" or "complex" PSM formulation
#' @return Vector of PPS survival function values
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
#'   )
#' # calc_surv_psmpps(totime=1:10,
#' #   fromtime=rep(1,10),
#' #   ptdata=bosonc,
#' #   dpam=params,
#' #   type="simple")
calc_surv_psmpps <- function(totime, fromtime=0, ptdata, dpam, psmtype="simple") {
  # Declare local variables
  cumH <- NULL
  # Hazard function
  hazfn <- function(x) {
    calc_haz_psm(timevar=x, dpam=dpam, ptdata=ptdata, psmtype=psmtype)$adj$pps
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

#' Obtain adjusted and unadjusted PSM hazards
#' @description EXPERIMENTAL. Obtain adjusted and unadjusted PSM hazards for given endpoint and time
#' @param endpoint Endpoint for which hazard is required (TTP, PPD, PFS, OS or PPS)
#' @inheritParams calc_haz_psm
#' @param psmtype Type of PSM - simple or complex
#' @importFrom rlang .data
#' @return `adj` is the hazard adjusted for constraints, `unadj` is the unadjusted hazard
pickout_psmhaz <- function(timevar, endpoint=NA, dpam, psmtype) {
  # Run calculation of all hazards
  allhaz <- calc_haz_psm(timevar, dpam, psmtype)
  # Required hazard, unadjusted
  h_unadj <- dplyr::case_when(
    endpoint=="TTP" ~ allhaz$unadj$ttp,
    endpoint=="PPD" ~ allhaz$unadj$ppd,
    endpoint=="PFS" ~ allhaz$unadj$pfs,
    endpoint=="OS" ~ allhaz$unadj$os,
    endpoint=="PPS" ~ allhaz$unadj$pps,
    .default = NA
  )
  # Required hazard, adjusted
  h_adj <- dplyr::case_when(
    endpoint=="TTP" ~ allhaz$adj$ttp,
    endpoint=="PPD" ~ allhaz$adj$ppd,
    endpoint=="PFS" ~ allhaz$adj$pfs,
    endpoint=="OS" ~ allhaz$adj$os,
    endpoint=="PPS" ~ allhaz$adj$pps,
    .default = NA
  )
  # Return adjusted and unadjusted hazard
  return(list(adj=h_adj, unadj=h_unadj, alladj=allhaz$adj))
}

#' Graph the PSM hazard functions
#' @description Graph the PSM hazard functions
#' @inheritParams pickout_psmhaz
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
#' @inherit pickout_psmhaz return
#' @importFrom rlang .data
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_par(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_par(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_par(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_par(fits$pfs, "aic")$fit,
#'   os = find_bestfit_par(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_par(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_par(fits$pps_cr, "aic")$fit
#' )
#' # Create graphics
#' # psmh_simple <- graph_psm_hazards(
#' #   timerange=(0:10)*6,
#' #   endpoint="OS",
#' #   ptdata=bosonc,
#' #   dpam=params,
#' #   psmtype="simple")
#' # psmh_simple$graph
graph_psm_hazards <- function(timevar, endpoint, ptdata, dpam, psmtype) {
  # Declare local variables
  Adjusted <- Unadjusted <- Time <- Hazard <- Method <- NULL
  # Convert endpoint to upper case text
  endpoint <- toupper(endpoint)
  # Pull out hazards to plot (inefficiently calls function twice, but is quite quick)
  adjhaz <- timevar |> purrr::map_vec(~pickout_psmhaz(.x, endpoint, dpam, psmtype)$adj)
  unadjhaz <- timevar |> purrr::map_vec(~pickout_psmhaz(.x, endpoint, dpam, psmtype)$unadj)
  # Create dataset for graphic
  result_data <- dplyr::tibble(Time=timevar, Adjusted=adjhaz, Unadjusted=unadjhaz) |>
      tidyr::pivot_longer(cols=c(Adjusted, Unadjusted),
                      names_to="Method", values_to="Hazard")
  # Draw graphic
  result_graph <- ggplot2::ggplot(result_data, ggplot2::aes(Time, Hazard)) +
    ggplot2::geom_line(ggplot2::aes(color = Method))
  # Return data and graphic
  return(list(data=result_data, graph=result_graph))
}

#' Graph the PSM survival functions
#' @description Graph the PSM survival functions
#' @inheritParams graph_psm_hazards
#' @inherit graph_psm_hazards return
#' @importFrom rlang .data
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_par(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit_par(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit_par(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit_par(fits$pfs, "aic")$fit,
#'   os = find_bestfit_par(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit_par(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit_par(fits$pps_cr, "aic")$fit
#' )
#' # Original OS graphic
#' graph_orig <- graph_survs(ptdata=bosonc, dpam=params)
#' graph_orig$graph$os
#' # New graphic illustrating effect of constraints on OS model
#' # psms_simple <- graph_psm_survs(
#' # timerange=6*(0:10),
#' # endpoint="OS",
#' # ptdata=bosonc,
#' # dpam=params,
#' # psmtype="simple")
#' # psms_simple$graph
graph_psm_survs <- function(timevar, endpoint, dpam, psmtype) {
  # Declare local variables
  Adjusted <- Unadjusted <- Time <- Survival <- Method <- NULL
  # Convert endpoint to upper case text
  endpoint <- toupper(endpoint)
  # Unadjusted hazard
  haz_unadj <- function(time) {
    pickout_psmhaz(time, endpoint, dpam, psmtype)$unadj
  }
  # Adjusted hazard
  haz_adj <- function(time) {
    pickout_psmhaz(time, endpoint, dpam, psmtype)$adj
  }
  # Unadjusted cumulative hazard
  cumhaz_unadj <- function(time) {
    intU <- stats::integrate(haz_unadj, lower=0, upper=time)
    if (intU$message=="OK") intU$value else NA
  }
  # Adjusted cumulative hazard
  cumhaz_adj <- function(time) {
    intA <- stats::integrate(haz_adj, lower=0, upper=time)
    if (intA$message=="OK") intA$value else NA
  }
  # Calculate cumulative hazards
  cumH_unadj <- seq(timevar) |> purrr::map_dbl(~cumhaz_unadj(.x))
  cumH_adj <- seq(timevar) |> purrr::map_dbl(~cumhaz_adj(.x))
  # Calculate survival values
  S_unadj <- exp(-cumH_unadj)
  S_adj <- exp(-cumH_adj)
  # Create dataset for graphic
  result_data <- dplyr::tibble(Time=timevar, Adjusted=S_adj, Unadjusted=S_unadj) |>
    tidyr::pivot_longer(cols=c(Adjusted, Unadjusted),
                        names_to="Method", values_to="Survival")
  # Draw graphic
  result_graph <- ggplot2::ggplot(result_data, ggplot2::aes(Time, Survival)) +
    ggplot2::geom_line(ggplot2::aes(color = Method))
  # Return data and graphic
  return(list(data=result_data, graph=result_graph))
}