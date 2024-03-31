# Discretized restricted mean durations


#' Restricted mean duration in progression-free for state transition models, discretized approximation
#' @description Calculates a discretized approximation for the mean duration in the progression-free state for both the state transition clock forward and clock reset models. Requires a carefully formatted list of fitted survival regressions for the necessary endpoints, and the time duration to calculate over.
#' This method is only valid for one-piece models, without lifetable adjustment.
#' @param dpam List of survival regressions for model endpoints. These must include time to progression (TTP) and pre-progression death (PPD).
#' @param Ty Time duration over which to calculate. Assumes input is in years, and patient-level data is recorded in weeks.
#' @param starting Vector of membership probabilities at time zero.
#' @param lifetable Optional. The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param discrate Discount rate (%) per year
#' @param state May be "PF", "PD" or "OS" being the progression-free or progressive disease states, or overall survival
#' @param model May be "PSM", "STM-CF" or "STM-CR" model structures
#' @return Numeric value in same time unit as patient-level data (weeks).
#' @include basics.R probgraphs.R
#' @seealso Full integral calculation in [prmd_pf_stm()]
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
#' discretized_rmd(dpam=params, state="pd", model="stm_cr")
#' rmd_pd_stm_cr(dpam=params)
discretized_rmd <- function(dpam, Ty=10, discrate=0, state, model, timestep=1) {
  # Declare local variables
  Tw <- tvec <- probs <- vn <- NULL
  # State and model should be lower case
  state <- tolower(state)
  model <- tolower(model)
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(1:floor(Tw/timestep)) + timestep/2
  # Vector of membership probabilities
  if (state=="pf") { 
    if (model=="psm") {probs <- prob_pf_psm(tvec, dpam)}
    else if (model=="stm_cf" | model=="stm_cr") {probs <- prob_pf_stm(tvec, dpam)}
  }
  else if (state=="pd") { 
    if (model=="psm") {probs <- prob_pd_psm(tvec, dpam)}
    else if (model=="stm_cf") {probs <- prob_pd_stm_cf(tvec, dpam)}
    else if (model=="stm_cr") {probs <- prob_pd_stm_cr(tvec, dpam)}
  }
  else if (state=="os") { 
    if (model=="psm") {probs <- prob_os_psm(tvec, dpam)}
    else if (model=="stm_cf") {probs <- prob_os_stm_cf(tvec, dpam)}
    else if (model=="stm_cr") {probs <- prob_os_stm_cr(tvec, dpam)}
  }
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Return value with starting adjustment
  sum(probs*vn) * timestep
}


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
#' drmd_psm(dpam=params, lifetable=lifetable)
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
  adjosprob <- constrain_survprob(osprob, lifetable=ltable, timevec=tvec)
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
#' drmd_stm_cf(dpam=params, lifetable=lifetable)
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
  adjsppd <- constrain_survprob(sppd, lifetable=ltable, timevec=tvec)
  adjos <- constrain_survprob(sos, lifetable=ltable, timevec=tvec)
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
#' drmd_stm_cr(dpam=params, lifetable=lifetable)
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
  adjsppd <- constrain_survprob(sppd, lifetable=ltable, timevec=tvec)
  adjos <- constrain_survprob(sos, lifetable=ltable, timevec=tvec)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(sttp*adjsppd*vn) * timestep
  os <- sum(adjos*vn) * timestep
  # Return values
  return(list(pf=pf, pd=os-pf, os=os))
}


