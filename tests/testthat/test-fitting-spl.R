# Testing file for fitting.R
# ==========================

# Fit_multi
# ---------

# Consider the two example datasets in the help files
# The bos dataset from the flexsurv package

bosonc <- create_dummydata("flexbosms")
mf1 <- fit_ends_mods_spl(bosonc)

# Fit "by hand"
knot_set <- c(1,1,1,2,2,2,3,3,3)
scale_set <- rep(c("hazard", "odds", "normal"),3)
Nmodels <- length(knot_set)

# PFS
ifits1_pfs <- purrr::map2(
  seq(knot_set), seq(scale_set),
    ~flexsurv::flexsurvspline(
      formula=survival::Surv(time=bosonc$pfs.durn,
                             event=bosonc$pfs.flag
      )~1,
      k = knot_set[.x],
      scale = scale_set[.y])
  )

# OS
ifits1_os <- purrr::map2(
  seq(knot_set), seq(scale_set),
  ~flexsurv::flexsurvspline(
    formula=survival::Surv(time=bosonc$os.durn,
                           event=bosonc$os.flag
    )~1,
    k = knot_set[.x],
    scale = scale_set[.y])
)

# TTP
ifits1_ttp <- purrr::map2(
  seq(knot_set), seq(scale_set),
  ~flexsurv::flexsurvspline(
    formula=survival::Surv(time=bosonc$ttp.durn,
                           event=bosonc$ttp.flag
    )~1,
    k = knot_set[.x],
    scale = scale_set[.y])
)

test_that("Fitted parameters, CI, SE, N and events match for bosonc data, TTP", {
  for (i in 1:Nmodels) {
    thismf <- mf1$ttp[[i]]
    thisif <- ifits1_ttp[[i]]
    expect_equal(thismf$result$aux$knots, thisif$aux$knots)
    expect_equal(thismf$result$aux$scale, thisif$aux$scale)
    expect_equal(thismf$result$N, thisif$N)
    expect_equal(thismf$result$events, thisif$events)
    expect_equal(thismf$result$res, thisif$res)
  }
})

test_that("Fitted parameters, CI, SE, N and events match for bosonc data, PFS", {
  for (i in 1:Nmodels) {
    thismf <- mf1$pfs[[i]]
    thisif <- ifits1_pfs[[i]]
    expect_equal(thismf$result$aux$knots, thisif$aux$knots)
    expect_equal(thismf$result$aux$scale, thisif$aux$scale)
    expect_equal(thismf$result$N, thisif$N)
    expect_equal(thismf$result$events, thisif$events)
    expect_equal(thismf$result$res, thisif$res)
  }
})

test_that("Fitted parameters, CI, SE, N and events match for bosonc data, OS", {
  for (i in 1:Nmodels) {
    thismf <- mf1$os[[i]]
    thisif <- ifits1_os[[i]]
    expect_equal(thismf$result$aux$knots, thisif$aux$knots)
    expect_equal(thismf$result$aux$scale, thisif$aux$scale)
    expect_equal(thismf$result$N, thisif$N)
    expect_equal(thismf$result$events, thisif$events)
    expect_equal(thismf$result$res, thisif$res)
  }
})

test_that("NA is produced when there are no distributions specified", {
  fitnull <- fit_ends_mods_spl(simdat=bosonc, k = NA, scale = NA)
  for (i in 1:6) {
    expect_equal(fitnull[[i]], NA)
  }
})

# findbest_survreg
# ----------------

# AIC

aics_ttp <- 1:9 |> purrr::map_vec(~mf1$ttp[[.x]]$result$AIC)
aics_pfs <- 1:9 |> purrr::map_vec(~mf1$pfs[[.x]]$result$AIC)
aics_os <- 1:9 |> purrr::map_vec(~mf1$os[[.x]]$result$AIC)
best_ttp <- which.min(aics_ttp)
best_pfs <- which.min(aics_pfs)
best_os <- which.min(aics_os)

test_that("findbest_survreg finds the best fits by min AIC", {
  expect_equal(find_bestfit_spl(mf1$ttp, "aic")$fit$AIC, aics_ttp[best_ttp])
  expect_equal(find_bestfit_spl(mf1$pfs, "aic")$fit$AIC, aics_pfs[best_pfs])
  expect_equal(find_bestfit_spl(mf1$os, "aic")$fit$AIC, aics_os[best_os])
  expect_equal(find_bestfit_spl(mf1$ttp, "aic")$fit$res[,1], mf1$ttp[[best_ttp]]$result$res[,1])
  expect_equal(find_bestfit_spl(mf1$pfs, "aic")$fit$res[,1], mf1$pfs[[best_pfs]]$result$res[,1])
  expect_equal(find_bestfit_spl(mf1$os, "aic")$fit$res[,1], mf1$os[[best_os]]$result$res[,1])
})

# BIC

ll_ttp <- 1:Nmodels |> purrr::map_vec(~mf1$ttp[[.x]]$result$loglik)
ll_pfs <- 1:Nmodels |> purrr::map_vec(~mf1$pfs[[.x]]$result$loglik)
ll_os <- 1:Nmodels |> purrr::map_vec(~mf1$os[[.x]]$result$loglik)
np_ttp <- 1:Nmodels |> purrr::map_vec(~mf1$ttp[[.x]]$result$npars)
np_pfs <- 1:Nmodels |> purrr::map_vec(~mf1$pfs[[.x]]$result$npars)
np_os <- 1:Nmodels |> purrr::map_vec(~mf1$os[[.x]]$result$npars)
N_ttp <- 1:Nmodels |> purrr::map_vec(~mf1$ttp[[.x]]$result$N)
N_pfs <- 1:Nmodels |> purrr::map_vec(~mf1$pfs[[.x]]$result$N)
N_os <- 1:Nmodels |> purrr::map_vec(~mf1$os[[.x]]$result$N)
bic_ttp <- np_ttp*log(N_ttp)-2*ll_ttp
bic_pfs <- np_pfs*log(N_pfs)-2*ll_pfs
bic_os <- np_os*log(N_os)-2*ll_os
bestb_ttp <- which.min(bic_ttp)
bestb_pfs <- which.min(bic_pfs)
bestb_os <- which.min(bic_os)

test_that("findbest_survreg finds the best fits by min BIC", {
  expect_equal(min(find_bestfit_spl(mf1$ttp, "bic")$results$bic), bic_ttp[bestb_ttp])
  expect_equal(min(find_bestfit_spl(mf1$pfs, "bic")$results$bic), bic_pfs[bestb_pfs])
  expect_equal(min(find_bestfit_spl(mf1$os, "bic")$results$bic), bic_os[bestb_os])
  expect_equal(find_bestfit_spl(mf1$ttp, "bic")$fit$res[,1], mf1$ttp[[bestb_ttp]]$result$res[,1])
  expect_equal(find_bestfit_spl(mf1$pfs, "bic")$fit$res[,1], mf1$pfs[[bestb_pfs]]$result$res[,1])
  expect_equal(find_bestfit_spl(mf1$os, "bic")$fit$res[,1], mf1$os[[bestb_os]]$result$res[,1])
})
