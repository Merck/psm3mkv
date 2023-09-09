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
# Functions relating to calculating Brier scores
# brier.R
# ==================================================================

#' Calculate Integrated Brier Scores for Overall Survival
#' @description Calculate integrated Brier score for overall survival for each decision model. This function is essentially a wrapper to [SurvMetrics::IBS] for use with this package.
#' @param ptdata Patient-level dataset
#' @param dpam List of survival regressions for model endpoints. Should include all six endpoints.
#' @return List containing
#' - osibs: Vector of integrated Brier scores for PSM, STM-CF and STM-CR decision models
#' - timerange: Vector of the time range over which the integrated Brier score is calculated
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
#' calc_ibs(bosonc, params)
calc_ibs <- function(ptdata, dpam) {
  cat("Initial calculations \n")
  # Pull out event times
  time_interest <- sort(ptdata$os.durn[ptdata$os.flag==1])
  time_interest <- time_interest[-which(time_interest==lag(time_interest))]
  # Pull out numbers of patients and event times
  n_pts <- length(ptdata$os.durn)
  n_times <- length(time_interest)
  # Basic survival object
  survobj <- sort(survival::Surv(ptdata$os.durn, ptdata$os.flag))
  # Matrices of predicted probabilities, by model
  cat("Calculating fitted probabilities \n")
  sp_mat_psm <- matrix(prob_os_psm(time_interest, dpam),
                        nrow=n_pts,
                        ncol=n_times,
                        byrow=TRUE)
  sp_mat_stm_cf <- matrix(prob_os_stm_cf(time_interest, dpam),
                       nrow=n_pts,
                       ncol=n_times,
                       byrow=TRUE)
  sp_mat_stm_cr <- matrix(prob_os_stm_cr(time_interest, dpam),
                       nrow=n_pts,
                       ncol=n_times,
                       byrow=TRUE)
  # IBS values
  cat("Calculating IBS values \n")
  ibs_psm <- SurvMetrics::IBS(survobj, sp_mat_psm, time_interest)
  ibs_stm_cf <- SurvMetrics::IBS(survobj, sp_mat_stm_cf, time_interest)
  ibs_stm_cr <- SurvMetrics::IBS(survobj, sp_mat_stm_cr, time_interest)
  # Pull together results
  osibs <- c(ibs_psm, ibs_stm_cf, ibs_stm_cr)
  names(osibs) <- c("PSM", "STM-CF", "STM-CR")
  timerange <- c(min(time_interest), max(time_interest))
  return(list(osibs=osibs, timerange=timerange))
}