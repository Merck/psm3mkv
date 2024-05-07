#' Compare likelihoods of PSMs
#'
#' @param ptdata 
#' @param fitslist 
#' @param cuttime 
#'
#' @return
#' @export
#'
#' @examples
#' # Fit parametric distributions to a dataset
#' bosonc <- create_dummydata("flexbosms")
#' parfits <- fit_ends_mods_par(bosonc)
#' # Present comparison of likelihood calculations
#' compare_psm_likes(bosonc, parfits)
compare_psm_likes <- function(ptdata, fitslist, cuttime=0) {
  
  # Number of distributions for each endpoint
  nttp <- length(fitslist$ttp)
  npfs <- length(fitslist$pfs)
  nos <- length(fitslist$os)
  
  # Best fits for each endpoint - AIC
  aic_best <- list(
    ttp = find_bestfit(fitslist$ttp, crit="aic")$fit$dlist$name,
    pfs = find_bestfit(fitslist$pfs, crit="aic")$fit$dlist$name,
    os = find_bestfit(fitslist$os, crit="aic")$fit$dlist$name
  )
  
  # Best fits for each endpoint - BIC
  bic_best <- list(
    ttp = find_bestfit(fitslist$ttp, crit="bic")$fit$dlist$name,
    pfs = find_bestfit(fitslist$pfs, crit="bic")$fit$dlist$name,
    os = find_bestfit(fitslist$os, crit="bic")$fit$dlist$name
  )
  
  res <- tibble::tibble(
    id = 1:(npfs*nos*(nttp+1)),
    ttp_meth = NA,
    pfs_dist = NA,
    os_dist = NA,
    ll = NA,
    npar = NA,
    npts = oslist[[1]]$result$N
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
  
  # PSM-simple
  message("PSM simple \n")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (p in 1:npfs) {
    thisfit$pfs <- fitslist$pfs[[p]]$result
    for (o in 1:nos) {
      thisfit$os <- fitslist$os[[o]]$result
      resrow <- (p-1)*nos + o
      res$ttp_meth[resrow] <- "simple"
      res$pfs_dist[resrow] <- thisfit$pfs$dlist$name
      res$os_dist[resrow] <- thisfit$os$dlist$name
      res$ll[resrow] <- slike_simple(thisfit)
      res$npar[resrow] <- thisfit$pfs$npars + thisfit$os$npars + 1
    }
  }
  
  # PSM-complex
  message("PSM complex \n")
  thisfit <- list(ttp=NA, pfs=NA, os=NA)
  for (t in 1:nttp) {
    thisfit$ttp <- fitslist$ttp[[t]]$result
    for (p in 1:npfs) {
      thisfit$pfs <- fitslist$pfs[[p]]$result
      for (o in 1:nos) {
        thisfit$os <- fitslist$os[[o]]$result
        resrow <- nos*npfs + (t-1)*npfs*nos + (p-1)*nos + o
        cat(".")
        res$ttp_meth[resrow] <- thisfit$ttp$dlist$name
        res$pfs_dist[resrow] <- thisfit$pfs$dlist$name
        res$os_dist[resrow] <- thisfit$os$dlist$name
        res$ll[resrow] <- slike_complex(thisfit)
        res$npar[resrow] <- thisfit$ttp$npars + thisfit$pfs$npars + thisfit$os$npars
      }
    }
  }
  
  # Add AIC and BIC, with ranks
  res <- res |>
    dplyr::mutate(
      aic = 2*npar-2*ll,
      bic = npar*log(npts)-2*ll,
      rank_aic = rank(aic),
      rank_bic = rank(bic)
    )
  
  return(list(results=res, aic=aic_best, bic=bic_best))
}