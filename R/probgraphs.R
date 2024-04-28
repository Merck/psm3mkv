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
# Functions that derive probabilities and create graphics
# probgraphs.R
# ==================================================================

# Probability functions:
#   prob_pf_psm: PF state, PSM model
#   prob_pf_stm: PF state, STM-CR or STM-CF models
#   prob_os_psm: OS states, PSM model
#   prob_pps: PD state based on time from progression (STM-CR model)
#   prob_pd_stm_cr: PD state, STM-CR model
#   prob_pd_stm_cf: PD state, STM-CF model
# Probability functions provided for completeness only:
#   prob_pd_psm: PD state, PSM model
#   prob_os_stm_cr: OS states, STM-CR model
#   prob_os_stm_cf: OS states, STM-CF model
# Graphing function: graph_survs

#' Calculate probability of being progression free in partitioned survival model
#' @description Calculates membership probability for the progression free state, at a particular time (vectorized), given a partitioned survival model with given statistical distributions and parameters.
#' @param time Time (numeric and vectorized)
#' @param dpam List of survival regressions for model endpoints. This must include progression-free survival (PFS).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pf_psm(0:100, params)
#' }
prob_pf_psm <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  survprob <- NULL
  # Calculate survival for each time
  survprob <- time |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$pfs))
  starting[1] * survprob
}

#' Calculate probability of being progression free in either state transition model (clock forward or clock reset)
#' @description Calculates membership probability for the progression free state, at a particular time (vectorized), given either state transition model (clock forward or clock reset) with given statistical distributions and parameters.
#' @param time Time (numeric and vectorized)
#' @param dpam List of survival regressions for model endpoints. This must include pre-progression death (PPD) and time to progression (TTP).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pf_stm(0:100, params)
#' }
prob_pf_stm <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  s1 <- s2 <- NULL
  # Calculate S_PPD and S_TTP
  s1 <- time |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ppd))
  s2 <- time |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ttp))
  starting[1] * s1 * s2
}

#' Calculate probability of being alive in a partitioned survival model
#' @description Calculates membership probability of being alive at a particular time (vectorized), given either state transition model (clock forward or clock reset) with given statistical distributions and parameters. This is the sum of membership probabilities in the progression free and progressed disease states.
#' @param time Time (numeric and vectorized)
#' @param dpam List of survival regressions for model endpoints. This must include overall survival (OS).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_os_psm(0:100, params)
#' }
prob_os_psm <- function(time, dpam, starting=c(1, 0, 0)){
  # Declare local variables
  survprob <- NULL
  # Calculate S_OS
  survprob <- time |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$os))
  (starting[1] + starting[2]) * survprob
}

#' Calculate probability of post progression survival under the state transition clock reset model
#' @description Calculates probability of post progression survival at a given time from progression (vectorized). This probability is from the state transition clock reset model, according to the given statistical distributions and parameters.
#' @param time Time (numeric and vectorized) from baseline - not time from progression.
#' @param dpam List of survival regressions for model endpoints. This must include post progression survival calculated under the clock reset state transition model.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pps_cr(0:100, params)
#' }
prob_pps_cr <- function(time, dpam) {
  time |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$pps_cr))
}

#' Calculate probability of post progression survival under the state transition clock forward model
#' @description Calculates probability of post progression survival at a given time from progression (vectorized). This probability is from the state transition clock forward model, according to the given statistical distributions and parameters.
#' @param ttptimes Time (numeric and vectorized) from progression - not time from baseline.
#' @param ppstimes Time (numeric and vectorized) of progression
#' @param dpam List of survival regressions for model endpoints. This must include post progression survival calculated under the clock forward state transition model.
#' @return Vector of the mean probabilities of post-progression survival at each PPS time, averaged over TTP times.
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pps_cf(0:100, 0:100, params)
#' }
prob_pps_cf <- function(ttptimes, ppstimes, dpam) {
  # Declare local variables
  s1 <- rel <- s2 <- meanrel <- durn <- NULL
  s1 <- ttptimes |>
    purrr::map_dbl(~calc_surv(.x, survobj=dpam$pps_cf))
  # Calculate expected survival for each patient, then average across patients
  rel <- s2 <-
    matrix(0, nrow=length(ppstimes), ncol=length(ttptimes))
  meanrel <- rep(0, length(ppstimes))
  for (i in seq_len(length(ppstimes))) {
    for(j in seq_len(length(ttptimes))) {
      durn <- ttptimes[j] + ppstimes[i]
      s2[i,j] <- calc_surv(durn, pps.ts$type, pps.ts$spec)
      rel[i,j] <- s2[i,j]/s1[j]
    }
    meanrel[i] <- mean(rel[i,1:length(ttptimes)])
  }
  meanrel
}

#' Calculate membership probability of progressed disease state in a partitioned survival model
#' @description Calculates membership probability of having progressed disease at a particular time (vectorized), given the partitioned survival model with certain statistical distributions and parameters.
#' @param time Time (numeric and vectorized)
#' @param dpam List of survival regressions for model endpoints. This must include progression-free survival (PFS) and overall survival (OS).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pd_psm(0:100, params)
#' }
prob_pd_psm <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  os <- pf <- NULL
  # Calculations
  os <- time |> purrr::map_dbl(~prob_os_psm(.x, dpam, starting))
  pf <- time |> purrr::map_dbl(~prob_pf_psm(.x, dpam, starting))
  os-pf
}

#' Calculate probability of having progressed disease under the state transition clock reset model
#' @description Calculates membership probability of the progressed disease state at a given time (vectorized). This probability is from the state transition clock reset model, according to the given statistical distributions and parameters.
#' @param time Time (numeric and vectorized) from baseline.
#' @param dpam List of survival regressions for model endpoints. This must include pre-progression death (PPD), time to progression (TTP) and post progression survival calculated under the clock reset model (PPS-CR).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pd_stm_cr(0:100, params)
#' }
prob_pd_stm_cr <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  ttp.ts <- ppd.ts <- pps.ts <- NULL
  int_pf <- int_pd <- NULL
  # Avoid integration if time==0
  if (time==0) {return(starting[2])}
  # Probability of PD, starting from PF
  integrand_pf <- function(u) {
    sttp <- calc_surv(u, survobj=dpam$ttp)
    sppd <- calc_surv(u, survobj=dpam$ppd)
    http <- calc_haz(u, survobj=dpam$ttp)
    spps <- calc_surv(time-u, survobj=dpam$pps_cr)
    sttp*sppd*http*spps
  }
  integrand_pf <- Vectorize(integrand_pf, "u")
  int_pf <- stats::integrate(integrand_pf, lower=0, upper=time)
  # Probability of PD, starting from PD
  int_pd <- calc_surv(time, pps.ts$type, pps.ts$spec)
  # Combined probability, given starting points
  starting[1] * int_pf$value + starting[2] * int_pd
}
prob_pd_stm_cr <- Vectorize(prob_pd_stm_cr, "time")

#' Calculate probability of having progressed disease under the state transition clock forward model
#' @description Calculates membership probability of the progressed disease state at a given time (vectorized). This probability is from the state transition clock forward model, according to the given statistical distributions and parameters.
#' @param time Time (numeric and vectorized) from baseline.
#' @param dpam List of survival regressions for model endpoints. This must include pre-progression death (PPD), time to progression (TTP) and post progression survival calculated under the clock forward model (PPS-CF).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_pd_stm_cf(0:100, params)
#' }
prob_pd_stm_cf <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  ttp.ts <- ppd.ts <- pps.ts <- NULL
  sppst <- int_pf <- int_pd <- NULL
  # Avoid integration if time==0
  if (time==0) {return(starting[2])}
  # SPPS
  sppst <- calc_surv(time, survobj=dpam$pps_cf)
  # Probability of PD, starting from PD
  integrand <- function(u) {
    sttp <- calc_surv(u, survobj=dpam$ttp)
    sppd <- calc_surv(u, survobj=dpam$ppd)
    http <- calc_haz(u, survobj=dpam$ttp)
    sppsu <- calc_surv(u, survobj=dpam$pps_cf)
    ifelse(sppsu==0, 0, sttp*sppd*http*sppst/sppsu)
  }
  integrand <- Vectorize(integrand, "u")
  int_pf <- stats::integrate(integrand, lower=0, upper=time)
  # Probability of PD, starting from PD
  int_pd <- sppst
  # Combined probability, given starting points
  starting[1] * int_pf$value + starting[2] * int_pd
}
prob_pd_stm_cf <- Vectorize(prob_pd_stm_cf, "time")

#' Calculate probability of being alive under the state transition clock reset model
#' @description Calculates membership probability of being alive at a given time (vectorized). This probability is from the state transition clock reset model, according to the given statistical distributions and parameters.
#' @param time Time (numeric and vectorized) from baseline.
#' @param dpam List of survival regressions for model endpoints. This must include pre-progression death (PPD), time to progression (TTP) and post progression survival calculated under the clock reset model (PPS-CR).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_os_stm_cr(0:100, params)
#' }
prob_os_stm_cr <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  pf <- pd <- NULL
  # Calculations
  pf <- time |> purrr::map_dbl(~prob_pf_stm(.x, dpam, starting))
  pd <- time |> purrr::map_dbl(~prob_pd_stm_cr(.x, dpam, starting))
  pf+pd
}

#' Calculate probability of being alive under the state transition clock forward model
#' @description Calculates membership probability of being alive at a given time (vectorized). This probability is from the state transition clock forward model, according to the given statistical distributions and parameters.
#' @param time Time (numeric and vectorized) from baseline.
#' @param dpam List of survival regressions for model endpoints. This must include pre-progression death (PPD), time to progression (TTP) and post progression survival calculated under the clock forward model (PPS-CF).
#' @param starting Vector of membership probabilities (PF, PD, death) at time zero.
#' @return Numeric value
#' @export
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_spl(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' prob_os_stm_cf(0:100, params)
#' }
prob_os_stm_cf <- function(time, dpam, starting=c(1, 0, 0)) {
  # Declare local variables
  pf <- pd <- NULL
  # Calculations
  pf <- time |> purrr::map_dbl(~prob_pf_stm(.x, dpam, starting))
  pd <- time |> purrr::map_dbl(~prob_pd_stm_cf(.x, dpam, starting))
  pf+pd
}

#' Graph the observed and fitted state membership probabilities
#' @description Graph the observed and fitted state membership probabilities for PF, PD, OS and PPS.
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
#' @param cuttime is the cut-off time for a two-piece model (default 0, indicating a one-piece model)
#' @return Four datasets and graphics as a list
#' @export
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_par(bosonc)
#' # Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#' )
#' # Create graphics
#' gs <- graph_survs(ptdata=bosonc, dpam=params)
#' gs$graph$pd
#' }
graph_survs <- function(ptdata, dpam, cuttime=0){
  cat("Creating KM \n")
  # Declare local variables
  ds <- dspps <- pfsfit <- osfit <- ppsfit <- NULL
  rnding <- kmpfs <- kmos <- kmpps <- kmdata <- gdata <- NULL
  time <- surv <- km_os <- km_pf <- endp <- timeplus <- NULL
  timecut <- cut_pf <- cut_os <- starting <- NULL
  # Calculations
  ds <- create_extrafields(ptdata, cuttime)
  dspps <- ds |> dplyr::filter(.data$pps.durn>0, .data$ttp.flag==1) 
  # Fit K-M before cutpoint
  pfsfit <- survival::survfit(
              survival::Surv(pfs.odurn, pfs.flag) ~ 1,
              data=ds)
  osfit <- survival::survfit(
              survival::Surv(os.odurn, os.flag) ~ 1,
              data=ds)
  ppsfit <- survival::survfit(
              survival::Surv(pps.odurn, pps.flag) ~ 1,
              data=dspps)
  # Obtain KM survivals at a single set of times
  kmtime <- sort(unique(c(pfsfit$time, osfit$time, ppsfit$time)))
  kmpfs <- summary(pfsfit, kmtime)
  kmos <- summary(osfit, kmtime)
  kmpps <- summary(ppsfit, kmtime)
  # Organize into a wide dataframe
  kmdata <- data.frame(
    time = c(kmpfs$time, kmos$time, kmpps$time),
    surv = c(kmpfs$surv, kmos$surv, kmpps$surv),
    epoint = c(rep("pf", length(kmpfs$time)), rep("os", length(kmos$time)), rep("pps", length(kmpps$time)))
  ) |>
    tidyr::pivot_wider(
      id_cols=time,
      names_prefix="km_",
      names_from="epoint",
      values_from=surv
      ) |>
    tibble::add_row(time=0, km_pf=1, km_os=1, km_pps=1) |>
    dplyr::mutate(km_pd = km_os-km_pf) |>
    dplyr::arrange(time)
  # Derive KM probabilities at cutpoint
  timecut <- max(kmdata$time[kmdata$time<=cuttime])
  cut_pf <- min(kmdata$km_pf[kmdata$time<=timecut])
  cut_os <- min(kmdata$km_os[kmdata$time<=timecut])
  starting <- c(cut_pf, cut_os-cut_pf, 1-cut_os)
  # Calculate fitted survival values
  cat("Calculating fitted curves \n")
  gdata <- kmdata |>
    dplyr::mutate(
      timeplus = pmax(0, .data$time-cuttime),
      psm_pf = prob_pf_psm(.data$timeplus, dpam, starting),
      psm_os = prob_os_psm(.data$timeplus, dpam, starting),
      psm_pd = .data$psm_os-.data$psm_pf,
      stmcf_pf = prob_pf_stm(.data$timeplus, dpam, starting),
      stmcr_pf = .data$stmcf_pf,
      stmcf_pd = prob_pd_stm_cf(.data$timeplus, dpam, starting),
      stmcr_pd = prob_pd_stm_cr(.data$timeplus, dpam, starting),
      stmcf_os = .data$stmcf_pf + .data$stmcf_pd,
      stmcr_os = .data$stmcr_pf + .data$stmcr_pd,
      stmcr_pps = prob_pps_cr(.data$time, dpam),
      stmcf_pps = prob_pps_cf(ttptimes=ds$ttp.durn, ppstimes=.data$time, dpam=dpam)
    ) |>
    dplyr::select(-timeplus) |>
    # Reshape into a long dataframe
    tidyr::pivot_longer(
      cols = !time,
      names_to = "surv"
    )
  # Pull out method and endpoint variables
  methep <- stringr::str_split(gdata$surv, "_", simplify=TRUE)
  gdata$method <- methep[,1]
  gdata$endp <- methep[,2]
  # Relabel stmcf and stmcr
  gdata$method[gdata$method=="stmcf"] <- "stm_cf"
  gdata$method[gdata$method=="stmcr"] <- "stm_cr"
  # Reshape wide by Method
  gdata <- gdata |>
    tidyr::pivot_wider(
      id_cols=c(time, endp),
      names_from="method",
      values_from="value"
    ) |>
    # Rename variables to: Time
    dplyr::rename(Time = time)
  # Set fitted values to NA before cuttime
  gdata$psm[gdata$Time < cuttime] <- NA
  gdata$stm_cr[gdata$Time < cuttime] <- NA
  gdata$stm_cf[gdata$Time < cuttime] <- NA
  # Internal function to draw graphic
  cat("Drawing plots \n")
  draw_2pgraphic <- function(graphds, xlabel="Time from baseline") {
    # Declare local variables
    endp <- Time <- Probability <- Method <- NULL
    # Reshape long by method
    longds <- graphds |>
      tidyr::pivot_longer(
                  cols=c("km", "psm", "stm_cf", "stm_cr"),
                  names_to="Method",
                  values_to="Probability"
                  ) |>
      dplyr::select(-endp)
    # Draw graphic
    ggplot2::ggplot(longds, ggplot2::aes(Time, Probability)) +
      ggplot2::geom_line(ggplot2::aes(color = Method)) +
      ggplot2::xlab(xlabel)
  }
  # Draw graphics
  graphlist <- list(
    pf = draw_2pgraphic(gdata[gdata$endp=="pf",]),
    pd = draw_2pgraphic(gdata[gdata$endp=="pd",]),
    os = draw_2pgraphic(gdata[gdata$endp=="os",]),
    pps = draw_2pgraphic(gdata[gdata$endp=="pps",], xlabel="Time from progression")
  )
  return(list(data=gdata, graph=graphlist))
}