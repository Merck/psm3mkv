# Testing for lhoods.R
# ====================

library(dplyr)
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
ll_psm <- calc_likes_psm(bosonc, params, cuttime=0)
llcomp_psm <- ll_psm$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_psm)
    )

ll_stmcf <- calc_likes_stm_cf(bosonc, params, cuttime=0)
llcomp_stmcf <- ll_stmcf$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_stm)
  )

ll_stmcr <- calc_likes_stm_cr(bosonc, params, cuttime=0)
llcomp_stmcr <- ll_stmcr$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_stm)
  )

ll_all <- calc_likes(bosonc, params, cuttime=0)

# Test likelihood values
test_that("Likelihood components match expected values - Par", {
  expect_equal(ll_psm$slikes$ll, llcomp_psm$slike)
  expect_equal(ll_stmcf$slikes$ll, llcomp_stmcf$slike)
  expect_equal(ll_stmcr$slikes$ll, llcomp_stmcr$slike)
})

test_that("Likelihood totals match expected values - Par", {
  expect_equal(ll_psm$ll, sum(llcomp_psm$slike))
  expect_equal(ll_stmcf$ll, sum(llcomp_stmcf$slike))
  expect_equal(ll_stmcr$ll, sum(llcomp_stmcr$slike))
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
ll_psm <- calc_likes_psm(bosonc, params_spl, cuttime=0)
llcomp_psm <- ll_psm$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_psm)
  )

ll_stmcf <- calc_likes_stm_cf(bosonc, params_spl, cuttime=0)
llcomp_stmcf <- ll_stmcf$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_stm)
  )

ll_stmcr <- calc_likes_stm_cr(bosonc, params_spl, cuttime=0)
llcomp_stmcr <- ll_stmcr$likedata |>
  group_by(outcome) |>
  summarize(
    slike=sum(llike_stm)
  )

ll_all <- calc_likes(bosonc, params_spl, cuttime=0)

# Test likelihood values
test_that("Likelihood components match expected values - Spl", {
  expect_equal(ll_psm$slikes$ll, llcomp_psm$slike)
  expect_equal(ll_stmcf$slikes$ll, llcomp_stmcf$slike)
  expect_equal(ll_stmcr$slikes$ll, llcomp_stmcr$slike)
})

test_that("Likelihood totals match expected values - Spl", {
  expect_equal(ll_psm$ll, sum(llcomp_psm$slike))
  expect_equal(ll_stmcf$ll, sum(llcomp_stmcf$slike))
  expect_equal(ll_stmcr$ll, sum(llcomp_stmcr$slike))
})
