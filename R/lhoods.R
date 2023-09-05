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
# Functions relating to calculating likelihoods
# lhoods.R
# ==================================================================

#' Obtain the type and specification as required in other package functions from a model fit
#' @param fitsurv Fitted model (either parametric or Royston-Parmar splines model)
#' @return List of type and specification
#' - `type` is "spl" for splines model or "par" for parametric model
#' - `spec` contains distribution (`dist`) and coefficients (`coefs`) if `type=="par"`
#' - `spec` contains gamma values (`gamma`), knot locations (log scale, `knots`) and scale (`scale`) for Royston-Parmar splines model, if `type=="spl"`
#' @examples
#' convert_fit2spec(dpam$pfs)
convert_fit2spec <- function(fitsurv) {
  par.dist <- fitsurv$dlist$name
  if (par.dist=="survspline") {
    type <- "spl"
    spl.gamma <- fitsurv$coefficients
    spl.knots <- fitsurv$aux$knots
    spl.scale <- fitsurv$aux$scale
    spec <- list(gamma=spl.gamma,
                 knots=spl.knots,
                 k=length(spl.knots)-2,
                 scale=spl.scale)
  } else {
    type <- "par"
    pars <- fitsurv$res[,1]
    spec <- list(dist=par.dist, pars=pars)
  }
  return(list(type=type, spec=spec))
}

#' Calculate likelihood for a simple three-state partitioned survival model
#' @description Calculate likelihood values and other summary output for a simple three-state partitioned survival model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].x three-state partitioned survival model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes]. Unlike [calc_likes_psm_complex], this likelihood function assumes a progression hazard can be derived from the PFS hazard function and the ratio of progression to PFS events from PF patients.
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
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
#' - `likedata`: Patient-level dataset with additional likelihood-related calculations.
#' - `coefsdists`: Summary table of distributions and parameters used for each endpoint.
#' - `slikes`: Total log-likelihood for each possible outcome
#' - `ll`: Total log-likelihood
#' - `params`: Number of parameters used in this model
#' - `AIC`: Akaike Information Criterion value for this model
#' - `BIC`: Bayesian Information Criterion value for this model
#' @seealso [calc_likes()], [calc_likes_psm_complex()], [calc_likes_stm_cf()], [calc_likes_stm_cr()]
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
#' calc_likes_psm_simple(bosonc, dpam=params)
calc_likes_psm_simple <- function(ptdata, dpam, cuttime=0) {
  # PFS
  pfs.ts <- convert_fit2spec(dpam$pfs)
  pfs.type <- pfs.ts$type
  pfs.spec <- pfs.ts$spec
  pfs.npars <- dpam$pfs$npars
  # OS
  os.ts <- convert_fit2spec(dpam$os)
  os.type <- os.ts$type
  os.spec <- os.ts$spec
  os.npars <- dpam$os$npars
  # Count parameters
  s_npars <- pfs.npars + os.npars + 1
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pmax(0, pfs.durn - cuttime),
      os.durn = pmax(0, os.durn - cuttime),
      ttp.durn = pmax(0, ttp.durn - cuttime)
    ) |>
    filter(pfs.durn>0) |>
    mutate(
      # Survival and hazard functions needed
      hpfsu = calc_haz(pfs.durn, pfs.type, pfs.spec),
      spfsu = calc_surv(pfs.durn, pfs.type, pfs.spec),
      hppdu = calc_haz_psm(timevar=pfs.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           type="simple")$pre,
      httpu = pmax(0, hpfsu-hppdu),
      sppstu = calc_surv_psmpps(totime=os.durn,
                                fromtime=pfs.durn,
                                ptdata=ptdata,
                                dpam=dpam,
                                type="simple"),
      hppst = calc_haz_psm(timevar=os.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           type="simple")$post
    )
  # Replace NA for zero hppst
  likedata$hppst[likedata$hppst==0] <- NA
  likedata <- likedata |>
    mutate(
      # Four possible outcomes
      f1 = (1-ttp.flag)*(1-os.flag),
      f2 = (1-ttp.flag)*os.flag,
      f3 = ttp.flag*(1-os.flag),
      f4 = ttp.flag*os.flag,
      # PSM likelihoods for each outcome
      # (u is PFS time rather than TTP time as in paper Table)
      # Add up and apply log
      like = case_when(
        f1==1 ~ spfsu,
        f2==1 ~ spfsu*hppdu,
        f3==1 ~ spfsu*httpu*sppstu,
        f4==1 ~ spfsu*httpu*sppstu*hppst,
        .default = NA
      ),
      llike = log(like),
      outcome = f1*1 + f2*2 + f3*3 + f4*4,
      chf = (f1+f2+f3+f4)==1,
      valid = is.na(llike)==FALSE
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
#' @export
#' @examples
#' bosonc <- create_dummydata("flexbosms")
#' fits <- fit_ends_mods_par(bosonc)
#' Pick out best distribution according to min AIC
#' params <- list(
#'   ppd = find_bestfit(fits$ppd, "aic")$fit,
#'   ttp = find_bestfit(fits$ttp, "aic")$fit,
#'   pfs = find_bestfit(fits$pfs, "aic")$fit,
#'   os = find_bestfit(fits$os, "aic")$fit,
#'   pps_cf = find_bestfit(fits$pps_cf, "aic")$fit,
#'   pps_cr = find_bestfit(fits$pps_cr, "aic")$fit
#'   )
#' calc_likes_psm_simple(bosonc, dpam=params)
calc_likes_psm_complex <- function(ptdata, dpam, cuttime=0) {
  # TTP
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ttp.type <- ttp.ts$type
  ttp.spec <- ttp.ts$spec
  ttp.npars <- dpam$ttp$npars
  # PFS
  pfs.ts <- convert_fit2spec(dpam$pfs)
  pfs.type <- pfs.ts$type
  pfs.spec <- pfs.ts$spec
  pfs.npars <- dpam$pfs$npars
  # OS
  os.ts <- convert_fit2spec(dpam$os)
  os.type <- os.ts$type
  os.spec <- os.ts$spec
  os.npars <- dpam$os$npars
  # Count parameters
  s_npars <- pfs.npars + os.npars + ttp.npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    dplyr::mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pmax(0, pfs.durn - cuttime),
      os.durn = pmax(0, os.durn - cuttime),
      ttp.durn = pmax(0, ttp.durn - cuttime)
    ) |>
    filter(pfs.durn>0) |>
    mutate(
      # Survival and hazard functions needed
      hpfsu = calc_haz(pfs.durn, pfs.type, pfs.spec),
      spfsu = calc_surv(pfs.durn, pfs.type, pfs.spec),
      hppdu = calc_haz_psm(timevar=pfs.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           type="complex")$pre,
      httpu = pmax(0, hpfsu-hppdu),
      sppstu = calc_surv_psmpps(totime=os.durn,
                                fromtime=pfs.durn,
                                ptdata=ptdata,
                                dpam=dpam,
                                type="complex"),
      hppst = calc_haz_psm(timevar=os.durn,
                           ptdata=ptdata,
                           dpam=dpam,
                           type="complex")$post
    )
  likedata$hppst[likedata$hppst==0] <- NA
  likedata <- likedata |>
    mutate(
      # Four possible outcomes
      f1 = (1-ttp.flag)*(1-os.flag),
      f2 = (1-ttp.flag)*os.flag,
      f3 = ttp.flag*(1-os.flag),
      f4 = ttp.flag*os.flag,
      # PSM likelihoods for each outcome
      # (u is PFS time rather than TTP time as in paper Table)
      # Add up and apply log
      like = case_when(
        f1==1 ~ spfsu,
        f2==1 ~ spfsu*hppdu,
        f3==1 ~ spfsu*httpu*sppstu,
        f4==1 ~ spfsu*httpu*sppstu*hppst,
        .default = NA
      ),
      llike = log(like),
      outcome = f1*1 + f2*2 + f3*3 + f4*4,
      chf = (f1+f2+f3+f4)==1,
      valid = is.na(llike)==FALSE
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

#' Calculate likelihood for a three-state clock forward markov model
#' @description Calculate likelihood values and other summary output for a three-state clock forward model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].
#' @inheritParams calc_likes_psm_simple
#' @inherit calc_likes_psm_simple return
#' @seealso [calc_likes()], [calc_likes_psm_simple()], [calc_likes_psm_complex()], [calc_likes_stm_cr()]
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
#' calc_likes_stm_cf(bosonc, dpam=params)
calc_likes_stm_cf <- function(ptdata, dpam, cuttime=0) {
  # Pull out distributions and parameters - PPD
  ppd.ts <- convert_fit2spec(dpam$ppd)
  ppd.type <- ppd.ts$type
  ppd.spec <- ppd.ts$spec
  ppd.npars <- dpam$ppd$npars
  # TTP
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ttp.type <- ttp.ts$type
  ttp.spec <- ttp.ts$spec
  ttp.npars <- dpam$ttp$npars
  # PPS_CF
  pps.ts <- convert_fit2spec(dpam$pps_cf)
  pps.type <- pps.ts$type
  pps.spec <- pps.ts$spec
  pps.npars <- dpam$pps_cf$npars
  # Count parameters
  s_npars <- ppd.npars + ttp.npars + pps.npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pfs.durn - cuttime,
      os.durn = os.durn - cuttime,
      ttp.durn = ttp.durn - cuttime
    ) |>
    filter(pfs.durn>0) |>
    mutate(
      # Survival and hazard functions needed
      httpu = calc_haz(pfs.durn, ttp.type, ttp.spec),
      sttpu = calc_surv(pfs.durn, ttp.type, ttp.spec),
      hppdu = calc_haz(pfs.durn, ppd.type, ppd.spec),
      sppdu = calc_surv(pfs.durn, ppd.type, ppd.spec),
      sppsu = calc_surv(pfs.durn, pps.type, pps.spec),
      sppst = calc_surv(os.durn, pps.type, pps.spec),
      hppst = calc_haz(os.durn, pps.type, pps.spec),
      # Four possible outcomes
      f1 = (1-ttp.flag)*(1-os.flag),
      f2 = (1-ttp.flag)*os.flag,
      f3 = ttp.flag*(1-os.flag),
      f4 = ttp.flag*os.flag,
      # STM likelihoods for each outcome
      # Add up and apply log
      like = case_when(
        f1==1 ~ sttpu*sppdu,
        f2==1 ~ sttpu*sppdu*hppdu,
        f3==1 ~ sttpu*sppdu*httpu*(sppst/sppsu),
        f4==1 ~ sttpu*sppdu*httpu*(sppst/sppsu)*hppst,
        .default = NA
      ),
      llike = log(like),
      outcome = f1*1 + f2*2 + f3*3 + f4*4,
      chf = (f1+f2+f3+f4)==1,
      valid = is.na(llike)==FALSE
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

#' Calculate likelihood for a three-state clock reset markov model
#' @description Calculate likelihood values and other summary output for a three-state clock reset model, given appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used). This function is called by [calc_likes].
#' @inheritParams calc_likes_psm_simple
#' @inherit calc_likes_psm_simple return
#' @seealso [calc_likes()], [calc_likes_stm_cf()], [calc_likes_psm_simple()], [calc_likes_psm_complex()]
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
#' calc_likes_stm_cr(bosonc, dpam=params)
calc_likes_stm_cr <- function(ptdata, dpam, cuttime=0) {
  # Pull out distributions and parameters
  # Pull out distributions and parameters - PPD
  ppd.ts <- convert_fit2spec(dpam$ppd)
  ppd.type <- ppd.ts$type
  ppd.spec <- ppd.ts$spec
  ppd.npars <- dpam$ppd$npars
  # TTP
  ttp.ts <- convert_fit2spec(dpam$ttp)
  ttp.type <- ttp.ts$type
  ttp.spec <- ttp.ts$spec
  ttp.npars <- dpam$ttp$npars
  # PPS_CR
  pps.ts <- convert_fit2spec(dpam$pps_cr)
  pps.type <- pps.ts$type
  pps.spec <- pps.ts$spec
  pps.npars <- dpam$pps_cr$npars
  # Count parameters
  s_npars <- ppd.npars + ttp.npars + pps.npars
  # Add on fields to dataset
  likedata <- tidyr::as_tibble(ptdata) |>
    mutate(
      # Adjust all durations for cutoff time
      pfs.durn = pfs.durn - cuttime,
      os.durn = os.durn - cuttime,
      ttp.durn = ttp.durn - cuttime
    ) |>
    filter(pfs.durn>0) |>
    mutate(
      # Survival and hazard functions needed
      httpu = calc_haz(pfs.durn, ttp.type, ttp.spec),
      sttpu = calc_surv(pfs.durn, ttp.type, ttp.spec),
      hppdu = calc_haz(pfs.durn, ppd.type, ppd.spec),
      sppdu = calc_surv(pfs.durn, ppd.type, ppd.spec),
      hppdt = calc_haz(os.durn, ppd.type, ppd.spec),
      sppstu = calc_surv(os.durn-pfs.durn, pps.type, pps.spec),
      hppstu = calc_haz(os.durn-pfs.durn, pps.type, pps.spec),
      # Four possible outcomes
      f1 = (1-ttp.flag)*(1-os.flag),   # No progression or death
      f2 = (1-ttp.flag)*os.flag,       # Death before progression
      f3 = ttp.flag*(1-os.flag),       # Progression, no death
      f4 = ttp.flag*os.flag,           # Progression, then death
      # STM likelihoods for each outcome
      # Add up and apply log
      like = case_when(
        f1==1 ~ sttpu*sppdu,
        f2==1 ~ sttpu*sppdu*hppdu,
        f3==1 ~ sttpu*sppdu*httpu*sppstu,
        f4==1 ~ sttpu*sppdu*httpu*sppstu*hppstu,
        .default = NA
      ),
      llike = log(like),
      outcome = f1*1 + f2*2 + f3*3 + f4*4,
      chf = (f1+f2+f3+f4)==1,
      valid = is.na(llike)==FALSE
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

#' Calculate likelihoods for three three-state models
#' @description Calculate likelihood values and other summary output for the following three state models: partitioned survival, clock forward markov, and clock reset markov. The function requires appropriately formatted patient-level data, a set of fitted survival regressions, and the time cut-off (if two-piece modeling is used).
#' @inheritParams calc_likes_psm_simple
#' @return Two outputs are returned:
#' 
#' `results` is a tibble of values and data relating to the likelihood for this model:
#' - `npts`: Number of patients analysed for each endpoint.
#' - `likedata`: Patient-level dataset with additional likelihood-related calculations.
#' - `coefsdists`: Summary table of distributions and parameters used for each endpoint.
#' - `slikes`: Total log-likelihood for each possible outcome
#' - `ll`: Total log-likelihood
#' - `params`: Number of parameters used in this model
#' - `AIC`: Akaike Information Criterion value for this model
#' - `BIC`: Bayesian Information Criterion value for this model
#' 
#' `llcomp` is a tibble providing a breakdown of the likelihood calculations by outcome. Outcomes are as follows:
#' - (1) refers to patients who remain alive and progression-free during the follow-up;
#' - (2) refers to patients who die without prior progression during the follow-up;
#' - (3) refers to patients who progress and then remain alive for the remaining follow-up, and
#' - (4) refers to patients who progress and die within the follow-up.
#'
#' The number of patients for each outcome are given for each model structure. You may confirm that these are identical across model structures.
#' The contribution of each patient group to the calculation of log-likelihood for each model is given in fields beginning `ll_`. This is helpful for understanding differences in likelihoods between model structures, according to patient outcomes.
#' @seealso This function calls [calc_likes_psm_simple()], [calc_likes_psm_complex()], [calc_likes_stm_cf()], and [calc_likes_stm_cr()].
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
#' calc_likes(bosonc, dpam=params)
calc_likes <- function(ptdata, dpam, cuttime=0) {
  methodnames <- c("psm_simple", "psm_complex", "stm_cf", "stm_cr")
  # Call PSM and STM functions
  list1 <- calc_likes_psm_simple(ptdata, dpam, cuttime)
  list2 <- calc_likes_psm_complex(ptdata, dpam, cuttime)
  list3 <- calc_likes_stm_cf(ptdata, dpam, cuttime)
  list4 <- calc_likes_stm_cr(ptdata, dpam, cuttime)
  # Create datasets for each method
  lldata1 <- list1$data |>
    select(ptid, outcome, llike, valid) |>
    mutate(methno=1)
  lldata2 <- list2$data |>
    select(ptid, outcome, llike, valid) |>
    mutate(methno=2)
  lldata3 <- list3$data |>
    select(ptid, outcome, llike, valid) |>
    mutate(methno=3)
  lldata4 <- list4$data |>
    select(ptid, outcome, llike, valid) |>
    mutate(methno=4)
  # Pull datasets together by columns
  llvdata <- lldata1 |>
    left_join(lldata2, by="ptid") |>
    left_join(lldata3, by="ptid") |>
    left_join(lldata4, by="ptid") |>
    mutate(validall = (valid.x*valid.y*valid.x.x*valid.y.y==1)) |>
    select(ptid, validall)
  # Pull datasets together by rows
  lldata <- lldata1 |>
    add_row(lldata2) |>
    add_row(lldata3) |>
    add_row(lldata4) |>
    mutate(methname = methodnames[methno]) |>
    left_join(llvdata, by="ptid")
  # Present likelihood by outcome, validity and model
  s1_long <- lldata |>
    summarise(
      npts = n(),
      ll = sum(llike),
      .by = c("valid", "outcome", "methno", "methname")
    )
  s1_wide <- s1_long |>
    select(-methname) |>
    tidyr::pivot_wider(names_from="methno",
                       names_prefix="ll_",
                       values_from="ll") |>
    arrange(valid, outcome)
  # Summarise across all patients by validity, irrespective of outcome
  s2_long <- lldata |>
    summarise(
      npts = n(),
      ll = sum(llike),
      .by = c("validall", "methno", "methname")
    )
  s2_wide <- s2_long |>
    tidyr::pivot_wider(names_from="validall",
                       values_from=c("npts", "ll"))
  # Summarise across all patients, irrespective of outcome or validity
  s3_long <- lldata |>
    summarise(
      npts = n(),
      ll = sum(llike),
      .by = c("validall", "methno", "methname")
    ) |>
  mutate(
    # Number of parameters
    nparam = case_when(
      methno == 1 ~ list1$npar,
      methno == 2 ~ list2$npar,
      methno == 3 ~ list3$npar,
      methno == 4 ~ list4$npar,
      .default = 0
    ),
    # Calculate AIC, BIC
    aic = ifelse(validall, 2*nparam-2*ll, NA),
    bic = ifelse(validall, nparam*log(npts)-2*ll, NA),
    # Calculate ranks, among pts where likelihood can be calculated
    rank_aic = rank(aic),
    rank_bic = rank(bic),
    rank_aic = ifelse(validall, rank_aic, NA),
    rank_bic = ifelse(validall, rank_bic, NA)
    ) |>
    arrange(validall, methno)
  # Return results
  return(list(
    detailed = s1_wide,
    valid = s2_wide,
    sumall = s3_long))
}