#' Compare likelihoods of PSMs
#' 
#' Compare the total log-likelihood values for the patient-level dataset after fitting PSM-simple and PSM-complex models to each combination of endpoint distributions
#' @inheritParams calc_allrmds
#' @param fitslist List of distribution fits to relevant endpoints, after calling `fit_ends_mods_par()` or `fit_ends_mods_spl()`
#' @importFrom rlang .data
#' @return List containing
#' - `res`: Dataset of calculation results for each model
#' - `ind_aic`: Set of statistical distributions for TTP, PFS and OS which individually fit each endpoint with the best (lowest) AIC
#' - `ind_bic`: Set of statistical distributions for TTP, PFS and OS which individually fit each endpoint with the best (lowest) BIC
#' - `jt_aic`:  Set of statistical distributions for TTP, PFS and OS which overall fit a PSM with the best (lowest) AIC
#' - `jt_bic`:  Set of statistical distributions for TTP, PFS and OS which overall fit a PSM with the best (lowest) BIC
#' @export
#' @examples
#' # Fit parametric distributions to a dataset
#' bosonc <- create_dummydata("flexbosms")
#' parfits <- fit_ends_mods_par(bosonc)
#' # Present comparison of likelihood calculations
#' compare_psm_likes(bosonc, parfits)
compare_psm_likes <- function(ptdata, fitslist, cuttime=0) {
  # Check that fitslist is a list of 6 endpoints
  if (length(fitslist)!=6) {stop("The list provided to fitslist must contain all 6 endpoints")}
  # Create local variables
  eps <- ndists <- aic_indbest <- bic_indbest <- bests <- res <- thisfit <- aic_jtbest <- bic_jtbest <- NULL
  # TTP, PFS and OS are endpoints 1, 3 and 4
  eps <- c(1, 3, 4)
  # Number of distributions for each endpoint
  ndists <- eps |>
    purrr::map_vec(~length(fitslist[[.x]]))
  # Best fits for each endpoint - AIC
  aic_indbest <- eps |>
    purrr::map_vec(~find_bestfit(fitslist[[.x]], crit="aic")$fit$dlist$name)
  # Best fits for each endpoint - BIC
  bic_indbest <- eps |>
    purrr::map_vec(~find_bestfit(fitslist[[.x]], crit="bic")$fit$dlist$name)
  # Join as a tibble
  bests <- rbind(aic_indbest, bic_indbest)
  bests <- tibble::tibble(
    ttp_meth = bests[,1],
    pfs_dist = bests[,2],
    os_dist = bests[,3],
    meth = "ind",
    ic = c("aic", "bic")
  )
  # Create results table for each model combination
  res <- tibble::tibble(
    id = 1:(ndists[3]*ndists[2]*(ndists[1]+1)),
    ttp_meth = NA,
    pfs_dist = NA,
    os_dist = NA,
    ll = NA,
    npar = NA,
    npts = fitslist$os[[1]]$result$N
  )
  # Create a safe calculation of the PSM-simple likelihood (returns NA on error)
  slike_simple <- purrr::possibly(
    ~calc_likes_psm_simple(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime)$ll[2],
    otherwise = NA)
  # Create a safe calculation of the PSM-complex likelihood (returns NA on error)
  slike_complex <- purrr::possibly(
    ~calc_likes_psm_complex(
      ptdata=ptdata,
      dpam=.x,
      cuttime=cuttime)$ll[2],
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
      res$ll[resrow] <- slike_simple(thisfit)
      res$npar[resrow] <- thisfit$pfs$npars + thisfit$os$npars + 1
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
        res$ll[resrow] <- slike_complex(thisfit)
        res$npar[resrow] <- thisfit$ttp$npars + thisfit$pfs$npars + thisfit$os$npars
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
  res$best_aic[res$ttp_meth==aic_indbest[1] & res$pfs_dist==aic_indbest[2] & res$os_dist==aic_indbest[3]] <- 1
  res$best_bic[res$ttp_meth==bic_indbest[1] & res$pfs_dist==bic_indbest[2] & res$os_dist==bic_indbest[3]] <- 1
  # Identify best distributions for overall AIC and BIC
  aic_jtbest <- res |>
    dplyr::filter(rank_aic==1) |>
    dplyr::select(ttp_meth, pfs_dist, os_dist) |>
    dplyr::mutate(meth="joint", ic="aic")
  bic_jtbest <- res |>
    dplyr::filter(rank_bic==1) |>
    dplyr::select(ttp_meth, pfs_dist, os_dist) |>
    dplyr::mutate(meth="joint", ic="aic")
  # Join together
  bests <- bests |>
    tibble::add_row(aic_jtbest) |>
    tibble::add_row(bic_jtbest)
  # Return
  return(list(results=res, bests=bests))
}