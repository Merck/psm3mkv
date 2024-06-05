#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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
# Functions relating to calculating likelihoods
# lhoods.R
# ==================================================================

#' Calculate likelihood for a simple three-state partitioned survival model
#' @description Calculate likelihood values and other summary output for a simple three-state partitioned survival model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].x three-state partitioned survival model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes]. Unlike [calc_likes_psm_complex], this likelihood function assumes a progression hazard can be derived from the PFS hazard function and the ratio of progression to PFS events from PF patients.
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - `ptid`: patient identifier
#' - `pfs.durn`: duration of PFS from baseline
#' - `pfs.flag`: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - `os.durn`: duration of OS from baseline
#' - `os.flag`: event flag for OS (=1 if death occurred, 0 for censoring)
#' - `ttp.durn`: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - `ttp.flag`: event flag for TTP (=1 if progression occurred, 0 for censoring).
#'
#' Survival data for all other endpoints (time to progression, pre-progression death, post-progression survival) are derived from PFS and OS.
#' @param dpam List of survival regressions for each endpoint:
#' - pre-progression death (PPD)
#' - time to progression (TTP)
#' - progression-free survival (PFS)
#' - overall survival (OS)
#' - post-progression survival clock forward (PPS-CF) and
#' - post-progression survival clock reset (PPS-CR).
#' @param cuttime Time cutoff - this is nonzero for two-piece models.
#' @return List of values and data relating to the likelihood for this model:
#' - `npts`: Number of patients analysed for each endpoint.
#' - `npar`: Number of parameters used in this model.
#' - `data`: A tibble of detailed likelihood calculations, where each row represents a patient.
#' @seealso [calc_likes()], [calc_likes_psm_complex()], [calc_likes_stm_cf()], [calc_likes_stm_cr()]
#' @importFrom rlang .data
#' @noRd
# Examples
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
#   )
# calc_likes_psm_simple(bosonc, dpam=params)
# }
calc_likes_psm_simple <- function(ptdata, dpam, cuttime=0) {
  # Declare local variablesL
  s_npars <- likedata <- npts <- llsum <- retlist <- pfs.durn <- NULL
  # Count parameters
  s_npars <- dpam$pfs$npars + dpam$os$npars + 1
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pmax(0, .data$pfs.durn - cuttime),
      os.durn = pmax(0, .data$os.durn - cuttime),
      ttp.durn = pmax(0, .data$ttp.durn - cuttime)
    ) |>
    dplyr::filter(pfs.durn>0) |>
    dplyr::mutate(
      # Survival and hazard functions needed
      hpfsu = calc_haz(.data$pfs.durn, survobj=dpam$pfs),
      spfsu = calc_surv(.data$pfs.durn, survobj=dpam$pfs),
      hppdu = calc_haz_psm(timevar=.data$pfs.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           psmtype="simple")$adj$ppd,
      httpu = pmax(0, .data$hpfsu-.data$hppdu),
      sppstu = calc_surv_psmpps(totime=.data$os.durn,
                                fromtime=.data$pfs.durn,
                                ptdata=ptdata,
                                dpam=dpam,
                                psmtype="simple"),
      hppst = calc_haz_psm(timevar=.data$os.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           psmtype="simple")$adj$pps
    )
  # Replace NA for zero hppst
  likedata$hppst[likedata$hppst==0] <- NA
  likedata <- likedata |>
    dplyr::mutate(
      # Four possible outcomes
      f1 = (1-.data$ttp.flag) * (1-.data$os.flag),
      f2 = (1-.data$ttp.flag) * .data$os.flag,
      f3 = .data$ttp.flag * (1-.data$os.flag),
      f4 = .data$ttp.flag * .data$os.flag,
      # PSM likelihoods for each outcome
      # (u is PFS time rather than TTP time as in paper Table)
      # Add up and apply log
      like = dplyr::case_when(
        f1==1 ~ spfsu,
        f2==1 ~ spfsu*hppdu,
        f3==1 ~ spfsu*httpu*sppstu,
        f4==1 ~ spfsu*httpu*sppstu*hppst,
        .default = NA
      ),
      llike = log(.data$like),
      outcome = .data$f1*1 + .data$f2*2 + .data$f3*3 + .data$f4*4,
      chf = (.data$f1+.data$f2+.data$f3+.data$f4)==1,
      valid = is.na(.data$llike)==FALSE
    )
  # Record other metrics
  npts <- c(length(likedata$ptid), length(likedata$ptid[likedata$valid==TRUE]))
  llsum <- c(sum(likedata$llike), sum(likedata$llike[likedata$valid==TRUE]))
  # Return list
  retlist <- list(npts = npts,
                  npar = s_npars,
                  ll = llsum,
                  data = likedata)
  return(retlist)
}

#' Calculate likelihood for a more complex three-state partitioned survival model
#' @description Calculate likelihood values and other summary output for a more complex three-state partitioned survival model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes()]. Unlike [calc_likes_psm_simple()], this likelihood function requires fitting to TTP.
#' @inheritParams calc_likes_psm_simple
#' @inherit calc_likes_psm_simple return
#' @seealso [calc_likes()], [calc_likes_psm_simple()], [calc_likes_psm_complex()], [calc_likes_stm_cr()]
#' @importFrom rlang .data
#' @noRd
# Examples
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_par(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_par(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_par(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_par(fits$pfs, "aic")$fit,
#   os = find_bestfit_par(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_par(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_par(fits$pps_cr, "aic")$fit
#   )
# calc_likes_psm_complex(bosonc, dpam=params)
calc_likes_psm_complex <- function(ptdata, dpam, cuttime=0) {
  # Declare local variables
  s_npars <- likedata <- npts <- llsum <- retlist <- pfs.durn <- NULL
  # Count parameters
  s_npars <- dpam$pfs$npars + dpam$os$npars + dpam$ttp$npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pmax(0, .data$pfs.durn - cuttime),
      os.durn = pmax(0, .data$os.durn - cuttime),
      ttp.durn = pmax(0, .data$ttp.durn - cuttime)
    ) |>
    dplyr::filter(pfs.durn>0) |>
    dplyr::mutate(
      # Survival and hazard functions needed
      hpfsu = calc_haz(.data$pfs.durn, survobj=dpam$pfs),
      spfsu = calc_surv(.data$pfs.durn, survobj=dpam$pfs),
      hppdu = calc_haz_psm(timevar=.data$pfs.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           psmtype="complex")$adj$ppd,
      httpu = pmax(0, .data$hpfsu-.data$hppdu),
      sppstu = calc_surv_psmpps(totime=.data$os.durn,
                                fromtime=.data$pfs.durn,
                                ptdata=ptdata,
                                dpam=dpam,
                                psmtype="complex"),
      hppst = calc_haz_psm(timevar=.data$os.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           psmtype="complex")$adj$pps
    )
  likedata$hppst[likedata$hppst==0] <- NA
  likedata <- likedata |>
    dplyr::mutate(
      # Four possible outcomes
      f1 = (1-.data$ttp.flag) * (1-.data$os.flag),
      f2 = (1-.data$ttp.flag) * .data$os.flag,
      f3 = .data$ttp.flag * (1-.data$os.flag),
      f4 = .data$ttp.flag * .data$os.flag,
      # PSM likelihoods for each outcome
      # (u is PFS time rather than TTP time as in paper Table)
      # Add up and apply log
      like = dplyr::case_when(
        f1==1 ~ spfsu,
        f2==1 ~ spfsu*hppdu,
        f3==1 ~ spfsu*httpu*sppstu,
        f4==1 ~ spfsu*httpu*sppstu*hppst,
        .default = NA
      ),
      llike = log(.data$like),
      outcome = .data$f1*1 + .data$f2*2 + .data$f3*3 + .data$f4*4,
      chf = (.data$f1+.data$f2+.data$f3+.data$f4)==1,
      valid = is.na(.data$llike)==FALSE
    )
  # Record other metrics
  npts <- c(length(likedata$ptid), length(likedata$ptid[likedata$valid==TRUE]))
  llsum <- c(sum(likedata$llike), sum(likedata$llike[likedata$valid==TRUE]))
  # Return list
  retlist <- list(npts = npts,
                  npar = s_npars,
                  ll = llsum,
                  data = likedata)
  return(retlist)
}

#' Calculate likelihood for a three-state clock forward state transition model
#' @description Calculate likelihood values and other summary output for a three-state clock forward state transition model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].
#' @inheritParams calc_likes_psm_simple
#' @inherit calc_likes_psm_simple return
#' @seealso [calc_likes()], [calc_likes_psm_simple()], [calc_likes_psm_complex()], [calc_likes_stm_cr()]
#' @importFrom rlang .data
#' @noRd
# Examples
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
#   )
# calc_likes_stm_cf(bosonc, dpam=params)
calc_likes_stm_cf <- function(ptdata, dpam, cuttime=0) {
  # Declare local variables
  s_npars <- likedata <- npts <- llsum <- retlist <- pfs.durn <- NULL
  # Count parameters
  s_npars <- dpam$ppd$npars + dpam$ttp$npars + dpam$pps_cf$npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = .data$pfs.durn - cuttime,
      os.durn = .data$os.durn - cuttime,
      ttp.durn = .data$ttp.durn - cuttime
    ) |>
    dplyr::filter(pfs.durn>0) |>
    dplyr::mutate(
      # Survival and hazard functions needed
      httpu = calc_haz(.data$pfs.durn, survobj=dpam$ttp),
      sttpu = calc_surv(.data$pfs.durn, survobj=dpam$ttp),
      hppdu = calc_haz(.data$pfs.durn, survobj=dpam$ppd),
      sppdu = calc_surv(.data$pfs.durn, survobj=dpam$ppd),
      sppsu = calc_surv(.data$pfs.durn, survobj=dpam$pps_cf),
      sppst = calc_surv(.data$os.durn, survobj=dpam$pps_cf),
      hppst = calc_haz(.data$os.durn, survobj=dpam$pps_cf),
      # Four possible outcomes
      f1 = (1-.data$ttp.flag) * (1-.data$os.flag),
      f2 = (1-.data$ttp.flag) * .data$os.flag,
      f3 = .data$ttp.flag * (1-.data$os.flag),
      f4 = .data$ttp.flag * .data$os.flag,
      # STM likelihoods for each outcome
      # Add up and apply log
      like = dplyr::case_when(
        f1==1 ~ sttpu*sppdu,
        f2==1 ~ sttpu*sppdu*hppdu,
        f3==1 ~ sttpu*sppdu*httpu*(sppst/sppsu),
        f4==1 ~ sttpu*sppdu*httpu*(sppst/sppsu)*hppst,
        .default = NA
      ),
      llike = log(.data$like),
      outcome = .data$f1*1 + .data$f2*2 + .data$f3*3 + .data$f4*4,
      chf = (.data$f1+.data$f2+.data$f3+.data$f4)==1,
      valid = is.na(.data$llike)==FALSE
    )
  # Pick out n and log like
  npts <- c(length(likedata$ptid),
            length(likedata$ptid[likedata$valid==TRUE]))
  llsum <- c(sum(likedata$llike),
              sum(likedata$llike[likedata$valid==TRUE]))
  # Return list
  retlist <- list(npts = npts,
                  npar = s_npars,
                  ll = llsum,
                  data = likedata)
  return(retlist)
}

#' Calculate likelihood for a three-state clock reset state transition model
#' @description Calculate likelihood values and other summary output for a three-state clock reset model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].
#' @inheritParams calc_likes_psm_simple
#' @inherit calc_likes_psm_simple return
#' @seealso [calc_likes()], [calc_likes_stm_cf()], [calc_likes_psm_simple()], [calc_likes_psm_complex()]
#' @importFrom rlang .data
#' @noRd
# Examples
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
#   )
# calc_likes_stm_cr(bosonc, dpam=params)
calc_likes_stm_cr <- function(ptdata, dpam, cuttime=0) {
  # Declare local variables
  s_npars <- likedata <- npts <- llsum <- retlist <- pfs.durn <- NULL
  # Count parameters
  s_npars <- dpam$ppd$npars + dpam$ttp$npars + dpam$pps_cr$npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = .data$pfs.durn - cuttime,
      os.durn = .data$os.durn - cuttime,
      ttp.durn = .data$ttp.durn - cuttime
    ) |>
    dplyr::filter(pfs.durn>0) |>
    dplyr::mutate(
      # Survival and hazard functions needed
      httpu = calc_haz(.data$pfs.durn, survobj=dpam$ttp),
      sttpu = calc_surv(.data$pfs.durn, survobj=dpam$ttp),
      hppdu = calc_haz(.data$pfs.durn, survobj=dpam$ppd),
      sppdu = calc_surv(.data$pfs.durn, survobj=dpam$ppd),
      hppdt = calc_haz(.data$os.durn, survobj=dpam$ppd),
      sppstu = calc_surv(.data$os.durn-pfs.durn, survobj=dpam$pps_cr),
      hppstu = calc_haz(.data$os.durn-.data$pfs.durn, survobj=dpam$pps_cr),
      # Four possible outcomes
      f1 = (1-.data$ttp.flag) * (1-.data$os.flag),   # No progression or death
      f2 = (1-.data$ttp.flag) * .data$os.flag,       # Death before progression
      f3 = .data$ttp.flag * (1-.data$os.flag),       # Progression, no death
      f4 = .data$ttp.flag * .data$os.flag,           # Progression, then death
      # STM likelihoods for each outcome
      # Add up and apply log
      like = dplyr::case_when(
        f1==1 ~ sttpu*sppdu,
        f2==1 ~ sttpu*sppdu*hppdu,
        f3==1 ~ sttpu*sppdu*httpu*sppstu,
        f4==1 ~ sttpu*sppdu*httpu*sppstu*hppstu,
        .default = NA
      ),
      llike = log(.data$like),
      outcome = .data$f1*1 + .data$f2*2 + .data$f3*3 + .data$f4*4,
      chf = (.data$f1+.data$f2+.data$f3+.data$f4)==1,
      valid = is.na(.data$llike)==FALSE
    )
  # Pick out n and log like
  npts <- c(length(likedata$ptid),
            length(likedata$ptid[likedata$valid==TRUE]))
  llsum <- c(sum(likedata$llike),
          sum(likedata$llike[likedata$valid==TRUE]))
  # Return list
  retlist <- list(npts = npts,
                  npar = s_npars,
                  ll = llsum,
                  data = likedata
                  )
  return(retlist)
}

#' Calculate likelihoods for three three-state model structures
#' @description Calculate likelihood values and other summary output for the following three state models structures: partitioned survival, clock forward state transition, and clock reset state transition. The function requires appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used).
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - `ptid`: patient identifier
#' - `pfs.durn`: duration of PFS from baseline
#' - `pfs.flag`: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - `os.durn`: duration of OS from baseline
#' - `os.flag`: event flag for OS (=1 if death occurred, 0 for censoring)
#' - `ttp.durn`: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - `ttp.flag`: event flag for TTP (=1 if progression occurred, 0 for censoring).
#'
#' Survival data for all other endpoints (time to progression, pre-progression death, post-progression survival) are derived from PFS and OS.
#' @param dpam List of survival regressions for each endpoint:
#' - pre-progression death (PPD)
#' - time to progression (TTP)
#' - progression-free survival (PFS)
#' - overall survival (OS)
#' - post-progression survival clock forward (PPS-CF) and
#' - post-progression survival clock reset (PPS-CR).
#' @param cuttime Time cutoff - this is nonzero for two-piece models.
#' @return A list of three tibbles:
#' `all` is a tibble of results for all patients:
#' - `methname`: the model structure or method.
#' - `npar`: is the number of parameters used by that method.
#' - `npts_1` to `npts_4` are the number of patients experiencing outcomes 1-4 respectively (see below), and `npts_tot` the total.
#' - `ll_1` to `ll_4` are the log-likelihood values for patients experiencing outcomes 1-4 respectively (see below), and `ll_tot` the total.
#' `valid` is a tibble of the same design as `all` but only in patients with valid likelihoods for all 4 methods
#' `sum` is a tibble in respect of patients with valid likelihoods for all 4 methods providing:
#' - `npts`: number of patients contributing results for this method.
#' - `npar`: number of parameters used by that method.
#' - `ll`: total log-likelihood
#' - `AIC`: Akaike Information Criterion value for this model
#' - `BIC`: Bayesian Information Criterion value for this model
#' 
#' The four outcomes are as follows:
#' - (1) refers to patients who remain alive and progression-free during the follow-up;
#' - (2) refers to patients who die without prior progression during the follow-up;
#' - (3) refers to patients who progress and then remain alive for the remaining follow-up, and
#' - (4) refers to patients who progress and die within the follow-up.
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
#'   )
#' calc_likes(bosonc, dpam=params)
#' }
calc_likes <- function(ptdata, dpam, cuttime=0) {
  # Declare local variables
  methodnames <- list1 <- list2 <- list3 <- list4 <- NULL
  lldata1 <- lldata2 <- lldata3 <- lldata4 <- methono <- methname <- NULL
  llvdata <- lldata <- npar <- npts_1 <- npts_2 <- npts_3 <- npts_4 <- npts_tot <- NULL
  ll_1 <- ll_2 <- ll_3 <- ll_4 <- ll_tot <- NULL
  ptid <- outcome <- llike <- valid <- valid.x <- valid.y <- valid.x.x <- valid.y.y <- validall <- NULL
  # methodnames <- c("psm_simple", "psm_complex", "stm_cf", "stm_cr")
  # Pull the likelihood calculations
  list1 <- calc_likes_psm_simple(ptdata, dpam, cuttime)
  list2 <- calc_likes_psm_complex(ptdata, dpam, cuttime)
  list3 <- calc_likes_stm_cf(ptdata, dpam, cuttime)
  list4 <- calc_likes_stm_cr(ptdata, dpam, cuttime)
  # Create datasets
  lldata1 <- list1$data |>
    dplyr::select(ptid, outcome, llike, valid) |>
    dplyr::mutate(
      methname = "psm_simple",
      npar = list1$npar
      )
  lldata2 <- list2$data |>
    dplyr::select(ptid, outcome, llike, valid) |>
    dplyr::mutate(
      methname = "psm_complex",
      npar = list2$npar
      )
  lldata3 <- list3$data |>
    dplyr::select(ptid, outcome, llike, valid) |>
    dplyr::mutate(
      methname = "stm_cf",
      npar = list3$npar
      )
  lldata4 <- list4$data |>
    dplyr::select(ptid, outcome, llike, valid) |>
    dplyr::mutate(
      methname = "stm_cr",
      npar = list4$npar
      )
  # Pull datasets together by columns
  llvdata <- lldata1 |>
    dplyr::left_join(lldata2, by="ptid") |>
    dplyr::left_join(lldata3, by="ptid") |>
    dplyr::left_join(lldata4, by="ptid") |>
    dplyr::mutate(validall = (valid.x*valid.y*valid.x.x*valid.y.y==1)) |>
    dplyr::select(ptid, validall)
  # Pull datasets together by rows
  lldata <- lldata1 |>
    dplyr::add_row(lldata2) |>
    dplyr::add_row(lldata3) |>
    dplyr::add_row(lldata4) |>
    dplyr::left_join(llvdata, by="ptid")
  # Detailed table, valid pts only
  tab_detval <- lldata |>
    dplyr::filter(validall==TRUE) |>
    dplyr::summarise(
      npts = dplyr::n(),
      npar = mean(npar),
      ll = sum(llike),
      .by = c("outcome", "methname")
    ) |>
    tidyr::pivot_wider(names_from="outcome",
                       values_from=c("npts", "ll")) |>
    dplyr::mutate(
      npts_tot = npts_1+npts_2+npts_3+npts_4,
      ll_tot = ll_1+ll_2+ll_3+ll_4,
    ) |>
    dplyr::select(methname, npar,
                   npts_1, npts_2, npts_3, npts_4, npts_tot,
                   ll_1, ll_2, ll_3, ll_4, ll_tot)
  # Detailed table, all patients (valid or not)
  tab_detall <- lldata |>
    dplyr::summarise(
      npts = dplyr::n(),
      npar = mean(npar),
      ll = sum(llike),
      .by = c("outcome", "methname")
    ) |>
    tidyr::pivot_wider(names_from="outcome",
                       values_from=c("npts", "ll")) |>
    dplyr::mutate(
      npts_tot = npts_1+npts_2+npts_3+npts_4,
      ll_tot = ll_1+ll_2+ll_3+ll_4,
    ) |>
    dplyr::select(methname, npar,
                  npts_1, npts_2, npts_3, npts_4, npts_tot,
                  ll_1, ll_2, ll_3, ll_4, ll_tot)
  # Totals, valid patients only
  tab_tots <- lldata |>
    dplyr::filter(validall==TRUE) |>
    dplyr::summarise(
      npts = dplyr::n(),
      npar = mean(npar),
      ll = sum(llike),
      .by = c("methname")
    ) |>
    dplyr::mutate(
      # Calculate AIC, BIC, ranks
      aic = 2*.data$npar - 2*.data$ll,
      bic = .data$npar * log(.data$npts) - 2 * .data$ll,
      rank_aic = rank(.data$aic),
      rank_bic = rank(.data$bic),
    )
    
  # Return results
  return(list(
    all = tab_detall,
    valid = tab_detval,
    sum = tab_tots))
}
