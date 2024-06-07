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
# Functions relating to comparing likelihoods of parametric and splines PSMs
# - compare_psm_likes()
# - compare_psm_likes_par()
# - compare_psm_likes_spl()
# ==================================================================

#' Compare likelihoods of PSMs
#' 
#' Compare the total log-likelihood values for the patient-level dataset after fitting PSM-simple and PSM-complex models to each combination of endpoint distributions
#' @inheritParams calc_allrmds
#' @param fitslist List of distribution fits to relevant endpoints, after calling `fit_ends_mods_par()` or `fit_ends_mods_spl()`
#' @importFrom rlang .data
#' @return List containing
#' - `results`: Dataset of calculation results for each model
#' - `bests`: Tibble indicating which is the best fitting model individually or jointly, to each endpoint, according to AIC or BIC
#' @export
#' @examples
#' # Fit parametric distributions to a dataset
#' bosonc <- create_dummydata("flexbosms")
#' parfits <- fit_ends_mods_par(bosonc)
#' \donttest{
#' splfits <- fit_ends_mods_spl(bosonc)
#' # Present comparison of likelihood calculations
#' compare_psm_likes(bosonc, parfits)
#' compare_psm_likes(bosonc, splfits)
#' }
compare_psm_likes <- function(ptdata, fitslist, cuttime=0) {
  # Determine whether fit is a spline, then call the right function
  if (fitslist$pfs[[1]]$result$dlist$name=="survspline") {
    compare_psm_likes_spl(ptdata, fitslist, cuttime)
  } else {
    compare_psm_likes_par(ptdata, fitslist, cuttime)
  }
}

# Compare likelihoods of PSMs for parametric models
compare_psm_likes_par <- function(ptdata, fitslist, cuttime=0) {
  # Check that fitslist is a list of 6 endpoints
  if (length(fitslist)!=6) {stop("The list provided to fitslist must contain all 6 endpoints")}
  # Create local variables
  eps <- ndists <- aic_indbest <- bic_indbest <- bests <- res <- thisfit <- aic_jtbest <- bic_jtbest <- NULL
  ll <- rank_aic <- ttp_meth <- pfs_dist <- os_dist <- rank_bic <- NULL
  eps <- c(1, 3, 4)
  # Best fits by AIC or BIC
  fits_aic <- eps |>
    purrr::map(~find_bestfit(fitslist[[.x]], crit="aic")$fit)
  fits_bic <- eps |>
    purrr::map(~find_bestfit(fitslist[[.x]], crit="bic")$fit)
  fits_all <- c(fits_aic, fits_bic)
  # Pull out names/distributions
  fits_names <- seq(2*length(eps)) |>
    purrr::map_vec(~fits_all[[.x]]$dlist$name)
  # Summarize parametric results in a table
  bests <- tibble::tibble(
    ttp_meth = fits_names[c(1, 4)],
    pfs_dist = fits_names[c(2, 5)],
    os_dist = fits_names[c(3, 6)],
    meth = "ind",
    ic = c("aic", "bic")
  )
  # Create results table for each model combination
  ndists <- eps |>
    purrr::map_vec(~length(fitslist[[.x]]))
  res <- tibble::tibble(
    id = 1:(ndists[3]*ndists[2]*(ndists[1]+1)),
    ttp_meth = NA,
    pfs_dist = NA,
    os_dist = NA,
    ll = NA,
    npar = NA,
    npts = NA
  )
  # Create a safe calculation of the PSM-simple likelihood (returns NA on error)
  slike_simple <- purrr::possibly(
    ~calc_likes_psm_simple(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime),
    otherwise = NA)
  # Create a safe calculation of the PSM-complex likelihood (returns NA on error)
  slike_complex <- purrr::possibly(
    ~calc_likes_psm_complex(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime),
    otherwise = NA)
  # Compute results for PSM-simple models
  message("Calculating PSM simple")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (p in 1:ndists[2]) {
    thisfit$pfs <- fitslist$pfs[[p]]$result
    for (o in 1:ndists[3]) {
      thisfit$os <- fitslist$os[[o]]$result
      resrow <- (p-1)*ndists[3] + o
      res$ttp_meth[resrow] <- "simple"
      res$pfs_dist[resrow] <- thisfit$pfs$dlist$name
      res$os_dist[resrow] <- thisfit$os$dlist$name
      psmsimple <- slike_simple(thisfit)
      if (is.na(psmsimple)[1]==FALSE) {
        res$ll[resrow] <- psmsimple$ll[2]
        res$npar[resrow] <- thisfit$pfs$npars + thisfit$os$npars + 1
        res$npts[resrow] <- psmsimple$npts[2]
      }
    }
  }
  # Compute results for PSM-complex models
  message("Calculating PSM complex")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (t in 1:ndists[1]) {
    thisfit$ttp <- fitslist$ttp[[t]]$result
    for (p in 1:ndists[2]) {
      thisfit$pfs <- fitslist$pfs[[p]]$result
      for (o in 1:ndists[3]) {
        thisfit$os <- fitslist$os[[o]]$result
        resrow <- t*ndists[3]*ndists[2] + (p-1)*ndists[3] + o
        res$ttp_meth[resrow] <- thisfit$ttp$dlist$name
        res$pfs_dist[resrow] <- thisfit$pfs$dlist$name
        res$os_dist[resrow] <- thisfit$os$dlist$name
        psmcomplex <- slike_complex(thisfit)
        if (is.na(psmcomplex)[1]==FALSE) {
          res$ll[resrow] <- psmcomplex$ll[2]
          res$npar[resrow] <- thisfit$ttp$npars + thisfit$pfs$npars + thisfit$os$npars
          res$npts[resrow] <- psmcomplex$npts[2]
        }
      }
    }
  }
  # Set log-likelihood values to NA if if cannot be calculated (=-Inf)
  res$ll[res$ll==-Inf] <- NA
  # Add AIC and BIC, with ranks
  message("Wrapping up")
  res <- res |>
    dplyr::mutate(
      aic = 2*.data$npar-2*ll,
      bic = .data$npar*log(.data$npts)-2*ll,
      rank_aic = rank(.data$aic),
      rank_bic = rank(.data$bic),
      best_aic = 0,
      best_bic = 0
    )
  # Identify best AIC and best BIC model
  res$best_aic[res$ttp_meth==bests$ttp_meth[1] & res$pfs_dist==bests$pfs_dist[1] & res$os_dist==bests$os_dist[1]] <- 1
  res$best_bic[res$ttp_meth==bests$ttp_meth[2] & res$pfs_dist==bests$pfs_dist[2] & res$os_dist==bests$os_dist[2]] <- 1
  # Identify best distributions for overall AIC and BIC
  aic_jtbest <- res |>
    dplyr::filter(rank_aic==1) |>
    dplyr::select(ttp_meth, pfs_dist, os_dist) |>
    dplyr::mutate(meth="joint", ic="aic")
  bic_jtbest <- res |>
    dplyr::filter(rank_bic==1) |>
    dplyr::select(ttp_meth, pfs_dist, os_dist) |>
    dplyr::mutate(meth="joint", ic="bic")
  # Join together
  bests <- bests |>
    tibble::add_row(aic_jtbest) |>
    tibble::add_row(bic_jtbest)
  # Return
  return(list(results=res, bests=bests))
}

# Compare likelihoods of PSMs for spline models
compare_psm_likes_spl <- function(ptdata, fitslist, cuttime=0) {
  # Check that fitslist is a list of 6 endpoints
  if (length(fitslist)!=6) {stop("The list provided to fitslist must contain all 6 endpoints")}
  # Create local variables
  eps <- ndists <- aic_indbest <- bic_indbest <- bests <- res <- thisfit <- aic_jtbest <- bic_jtbest <- NULL
  ll <- rank_aic <- ttp_meth <- ttp_knots <- pfs_scales <- pfs_knots <- os_scales <- os_knots <- rank_bic <- NULL
  # Best fits by AIC or BIC
  eps <- c(1, 3, 4)
  fits_aic <- eps |>
    purrr::map(~find_bestfit(fitslist[[.x]], crit="aic")$fit)
  fits_bic <- eps |>
    purrr::map(~find_bestfit(fitslist[[.x]], crit="bic")$fit)
  fits_all <- c(fits_aic, fits_bic)
  # Pull out scales and knots
  fits_scales <- seq(2*length(eps)) |>
    purrr::map_vec(~fits_all[[.x]]$aux$scale)
  fits_knots <- seq(2*length(eps)) |>
    purrr::map_vec(~length(fits_all[[.x]]$aux$knots))
  # Summarize parametric results in a table
  bests <- tibble::tibble(
    ttp_meth = fits_scales[c(1, 4)],
    ttp_knots = fits_knots[c(1, 4)],
    pfs_scales = fits_scales[c(2, 5)],
    pfs_knots = fits_knots[c(2, 5)],
    os_scales = fits_scales[c(3, 6)],
    os_knots = fits_knots[c(3, 6)],
    meth = "ind",
    ic = c("aic", "bic")
  )
  # Create results table for each model combination
  ndists <- eps |>
    purrr::map_vec(~length(fitslist[[.x]]))
  res <- tibble::tibble(
    id = 1:(ndists[3]*ndists[2]*(ndists[1]+1)),
    ttp_meth = NA,
    ttp_knots = NA,
    pfs_scales = NA,
    pfs_knots = NA,
    os_scales = NA,
    os_knots = NA,
    ll = NA,
    npar = NA,
    npts = NA
  )
  # Create a safe calculation of the PSM-simple likelihood (returns NA on error)
  slike_simple <- purrr::possibly(
    ~calc_likes_psm_simple(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime),
    otherwise = NA)
  # Create a safe calculation of the PSM-complex likelihood (returns NA on error)
  slike_complex <- purrr::possibly(
    ~calc_likes_psm_complex(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime),
    otherwise = NA)
  # Compute results for PSM-simple models
  message("Calculating PSM simple")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (p in 1:ndists[2]) {
    thisfit$pfs <- fitslist$pfs[[p]]$result
    for (o in 1:ndists[3]) {
      thisfit$os <- fitslist$os[[o]]$result
      resrow <- (p-1)*ndists[3] + o
      res$ttp_meth[resrow] <- "simple"
      res$ttp_knots[resrow] <- 0
      res$pfs_scales[resrow] <- thisfit$pfs$aux$scale
      res$pfs_knots[resrow] <- length(thisfit$pfs$aux$knots)
      res$os_scales[resrow] <- thisfit$os$aux$scale
      res$os_knots[resrow] <- length(thisfit$os$aux$knots)
      psmsimple <- slike_simple(thisfit)
      if (is.na(psmsimple)[1]==FALSE) {
        res$ll[resrow] <- psmsimple$ll[2]
        res$npar[resrow] <- thisfit$pfs$npars + thisfit$os$npars + 1
        res$npts[resrow] <- psmsimple$npts[2]
      }
    }
  }
  # Compute results for PSM-complex models
  message("Calculating PSM complex")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (t in 1:ndists[1]) {
    thisfit$ttp <- fitslist$ttp[[t]]$result
    for (p in 1:ndists[2]) {
      thisfit$pfs <- fitslist$pfs[[p]]$result
      for (o in 1:ndists[3]) {
        thisfit$os <- fitslist$os[[o]]$result
        resrow <- t*ndists[3]*ndists[2] + (p-1)*ndists[3] + o
        res$ttp_meth[resrow] <- thisfit$ttp$aux$scale
        res$ttp_knots[resrow] <- length(thisfit$ttp$aux$knots)
        res$pfs_scales[resrow] <- thisfit$pfs$aux$scale
        res$pfs_knots[resrow] <- length(thisfit$pfs$aux$knots)
        res$os_scales[resrow] <- thisfit$os$aux$scale
        res$os_knots[resrow] <- length(thisfit$os$aux$knots)
        psmcomplex <- slike_complex(thisfit)
        if (is.na(psmcomplex)[1]==FALSE) {
          res$ll[resrow] <- psmcomplex$ll[2]
          res$npar[resrow] <- thisfit$ttp$npars + thisfit$pfs$npars + thisfit$os$npars
          res$npts[resrow] <- psmcomplex$npts[2]
        }
      }
    }
  }
  # Set log-likelihood values to NA if if cannot be calculated (=-Inf)
  res$ll[res$ll==-Inf] <- NA
  # Add AIC and BIC, with ranks
  message("Wrapping up")
  res <- res |>
    dplyr::mutate(
      aic = 2*.data$npar-2*ll,
      bic = .data$npar*log(.data$npts)-2*ll,
      rank_aic = rank(.data$aic),
      rank_bic = rank(.data$bic),
      best_aic = 0,
      best_bic = 0
    )
  # Identify best AIC and best BIC model
  res$best_aic[res$ttp_meth==bests$ttp_meth[1] & res$ttp_knots==bests$ttp_knots[1] & res$pfs_scales==bests$pfs_scales[1] & res$pfs_knots==bests$pfs_knots[1] & res$os_scales==bests$os_scales[1] & res$os_knots==bests$os_knots[1]] <- 1
  res$best_bic[res$ttp_meth==bests$ttp_meth[2] & res$ttp_knots==bests$ttp_knots[2] & res$pfs_scales==bests$pfs_scales[2] & res$pfs_knots==bests$pfs_knots[2] & res$os_scales==bests$os_scales[2] & res$os_knots==bests$os_knots[2]] <- 1
  # Identify best distributions for overall AIC and BIC
  aic_jtbest <- res |>
    dplyr::filter(rank_aic==1) |>
    dplyr::select(ttp_meth, ttp_knots, pfs_scales, pfs_knots, os_scales, os_knots) |>
    dplyr::mutate(meth="joint", ic="aic")
  bic_jtbest <- res |>
    dplyr::filter(rank_bic==1) |>
    dplyr::select(ttp_meth, ttp_knots, pfs_scales, pfs_knots, os_scales, os_knots) |>
    dplyr::mutate(meth="joint", ic="bic")
  # Join together
  bests <- bests |>
    tibble::add_row(aic_jtbest) |>
    tibble::add_row(bic_jtbest)
  # Return
  return(list(results=res, bests=bests))
}