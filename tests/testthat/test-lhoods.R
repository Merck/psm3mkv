# Testing for lhoods.R
# ====================

bosonc <- create_dummydata("flexbosms")

# Parametric
# ----------

alldists <- c("exp", "weibullPH", "weibull", "llogis", "lnorm", "gamma",
              "gompertz","gengamma","gengamma.orig")

# Fit all distributions to all endpoints (except gen gamma for PPD)
af1 <- fit_ends_mods_par(bosonc,
                        cuttime = 0,
                        ppd.dist = alldists[1:7],
                        ttp.dist = alldists,
                        pfs.dist = alldists,
                        os.dist = alldists,
                        pps_cr.dist = alldists,
                        pps_cf.dist = alldists)

# Pick out the best fit for each endpoint
fit.ppd <- find_bestfit_par(af1$ppd, "aic")
fit.ttp <- find_bestfit_par(af1$ttp, "aic")
fit.pfs <- find_bestfit_par(af1$pfs, "aic")
fit.os <- find_bestfit_par(af1$os, "aic")
fit.pps_cf <- find_bestfit_par(af1$pps_cf, "aic")
fit.pps_cr <- find_bestfit_par(af1$pps_cr, "aic")

# Bring together
params <- list(ppd=fit.ppd$fit,
               ttp=fit.ttp$fit,
               pfs=fit.pfs$fit,
               os=fit.os$fit,
               pps_cf=fit.pps_cf$fit,
               pps_cr=fit.pps_cr$fit)

# Calculate likelihoods

# PSM simple
ll_psm_simple <- calc_likes_psm_simple(bosonc, params, cuttime=0)
llcomp_psm_simple <- ll_psm_simple$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
    )

# PSM complex
ll_psm_complex <- calc_likes_psm_complex(bosonc, params, cuttime=0)
llcomp_psm_complex <- ll_psm_complex$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

# STM CF
ll_stmcf <- calc_likes_stm_cf(bosonc, params, cuttime=0)
llcomp_stmcf <- ll_stmcf$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

# STM CR
ll_stmcr <- calc_likes_stm_cr(bosonc, params, cuttime=0)
llcomp_stmcr <- ll_stmcr$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

ll_all <- calc_likes(bosonc, params, cuttime=0)

# Test likelihood values
test_that("Likelihood totals match for parametric", {
  expect_equal(
    ll_psm_simple$ll[2],
    sum(llcomp_psm_simple$slike[llcomp_psm_simple$valid==TRUE])
    )
  expect_equal(
    ll_psm_complex$ll[2],
    sum(llcomp_psm_complex$slike[llcomp_psm_complex$valid==TRUE])
    )
  expect_equal(
    ll_stmcf$ll[2],
    sum(llcomp_stmcf$slike[llcomp_stmcf$valid==TRUE])
    )
  expect_equal(
    ll_stmcr$ll[2],
    sum(llcomp_stmcr$slike[llcomp_stmcr$valid==TRUE])
    )
})

# Splines
# ----------

# Fit all distributions to all endpoints (except gen gamma for PPD)
af2 <- fit_ends_mods_spl(bosonc)

# Pick out the best fit for each endpoint
fit.ppd <- find_bestfit_par(af2$ppd, "aic")
fit.ttp <- find_bestfit_par(af2$ttp, "aic")
fit.pfs <- find_bestfit_par(af2$pfs, "aic")
fit.os <- find_bestfit_par(af2$os, "aic")
fit.pps_cf <- find_bestfit_par(af2$pps_cf, "aic")
fit.pps_cr <- find_bestfit_par(af2$pps_cr, "aic")

# Bring together
params_spl <- list(ppd=fit.ppd$fit,
               ttp=fit.ttp$fit,
               pfs=fit.pfs$fit,
               os=fit.os$fit,
               pps_cf=fit.pps_cf$fit,
               pps_cr=fit.pps_cr$fit)

# Calculate likelihoods
# PSM simple
ll_psm_simple <- calc_likes_psm_simple(bosonc, params_spl, cuttime=0)
llcomp_psm_simple <- ll_psm_simple$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

# PSM complex
ll_psm_complex <- calc_likes_psm_complex(bosonc, params_spl, cuttime=0)
llcomp_psm_complex <- ll_psm_complex$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

# STM CF
ll_stmcf <- calc_likes_stm_cf(bosonc, params_spl, cuttime=0)
llcomp_stmcf <- ll_stmcf$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

# STM CR
ll_stmcr <- calc_likes_stm_cr(bosonc, params_spl, cuttime=0)
llcomp_stmcr <- ll_stmcr$data |>
  dplyr::group_by(valid, outcome) |>
  dplyr::summarize(
    slike=sum(llike)
  )

ll_all <- calc_likes(bosonc, params_spl, cuttime=0)

# Test likelihood values
test_that("Likelihood totals match for splines", {
  expect_equal(
    ll_psm_simple$ll[2],
    sum(llcomp_psm_simple$slike[llcomp_psm_simple$valid==TRUE])
  )
  expect_equal(
    ll_psm_complex$ll[2],
    sum(llcomp_psm_complex$slike[llcomp_psm_complex$valid==TRUE])
  )
  expect_equal(
    ll_stmcf$ll[2],
    sum(llcomp_stmcf$slike[llcomp_stmcf$valid==TRUE])
  )
  expect_equal(
    ll_stmcr$ll[2],
    sum(llcomp_stmcr$slike[llcomp_stmcr$valid==TRUE])
  )
})
