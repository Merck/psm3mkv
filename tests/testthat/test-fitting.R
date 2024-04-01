# Testing file for fitting.R
# ==========================

# Fit_multi
# ---------

# Consider the two example datasets in the help files
# The bos dataset from the flexsurv package

bosonc <- create_dummydata("flexbosms")

# Consider the two example datasets in the help files
# The bos dataset from the flexsurv package
alldists <- c("exp", "weibullPH", "weibull", "llogis", "lnorm", "gamma",
              "gompertz","gengamma")
mf1 <- fit_ends_mods_par(bosonc,
                         ppd.dist=alldists,
                         ttp.dist=alldists,
                         pfs.dist=alldists,
                         os.dist=alldists,
                         pps_cf.dist=alldists,
                         pps_cr.dist=alldists)

# Fit "by hand"
ifits_pfs <- alldists |>
  purrr::map(
    ~flexsurv::flexsurvreg(
      formula=survival::Surv(time=bosonc$pfs.durn,
                 event=bosonc$pfs.flag
                 )~1,
    dist=.x)
    )

ifits_os <- alldists |>
  purrr::map(
    ~flexsurv::flexsurvreg(
      formula=survival::Surv(time=bosonc$os.durn,
                             event=bosonc$os.flag
      )~1,
      dist=.x)
  )

ifits_ttp <- alldists |>
  purrr::map(
    ~flexsurv::flexsurvreg(
      formula=survival::Surv(time=bosonc$ttp.durn,
                             event=bosonc$ttp.flag
      )~1,
      dist=.x)
  )

test_that("Fitted parameters, CI, SE, N and events match for TTP", {
  for (i in seq(alldists)) {
    thismf <- mf1$ttp[[i]]$result
    thisif <- ifits_ttp[[i]]
    expect_equal(thismf$N, thisif$N)
    expect_equal(thismf$events, thisif$events)
    expect_equal(thismf$res, thisif$res, tolerance=0.001)
  }
})

test_that("Fitted parameters, CI, SE, N and events match for PFS", {
  for (i in seq(alldists)) {
    thismf <- mf1$pfs[[i]]$result
    thisif <- ifits_pfs[[i]]
    expect_equal(thismf$N, thisif$N)
    expect_equal(thismf$events, thisif$events)
    expect_equal(thismf$res, thisif$res, tolerance=0.001)
  }
})

test_that("Fitted parameters, CI, SE, N and events match for OS", {
  for (i in seq(alldists)) {
    thismf <- mf1$os[[i]]$result
    thisif <- ifits_os[[i]]
    expect_equal(thismf$N, thisif$N)
    expect_equal(thismf$events, thisif$events)
    expect_equal(thismf$res, thisif$res, tolerance=0.001)
  }
})

fitnull <- fit_ends_mods_par(simdat=bosonc,
                             ppd.dist=NA,
                             ttp.dist=NA,
                             pfs.dist=NA,
                             os.dist=NA,
                             pps_cf.dist=NA,
                             pps_cr.dist=NA)

test_that("NA is produced when there are no distributions specified", {
  for (i in 1:6) {
    expect_equal(fitnull[[i]], NA)
  }
})

# Check_posdef
# ------------

test_that("NA is produced when there are no distributions specified", {
  for (i in seq(alldists)) {
      expect_equal(
        det(chol(mf1$ttp[[i]]$result$opt$hessian))>0,
        check_posdef(mf1$ttp[[i]]$result)
      )
    expect_equal(
      det(chol(mf1$pfs[[i]]$result$opt$hessian))>0,
      check_posdef(mf1$pfs[[i]]$result)
    )
    expect_equal(
      det(chol(mf1$os[[i]]$result$opt$hessian))>0,
      check_posdef(mf1$os[[i]]$result)
    )
  }
})


# findbest_survreg
# ----------------

Ndists <- length(alldists)
aics_ttp <- 1:Ndists |> purrr::map_vec(~mf1$ttp[[.x]]$result$AIC)
aics_pfs <- 1:Ndists |> purrr::map_vec(~mf1$pfs[[.x]]$result$AIC)
aics_os <- 1:Ndists |> purrr::map_vec(~mf1$os[[.x]]$result$AIC)
best_ttp <- which.min(aics_ttp)
best_pfs <- which.min(aics_pfs)
best_os <- which.min(aics_os)

test_that("findbest_survreg finds the best fits by min AIC", {
  expect_equal(find_bestfit_par(mf1$ttp, "aic")$fit$AIC, aics_ttp[best_ttp])
  expect_equal(find_bestfit_par(mf1$pfs, "aic")$fit$AIC, aics_pfs[best_pfs])
  expect_equal(find_bestfit_par(mf1$os, "aic")$fit$AIC, aics_os[best_os])
  expect_equal(find_bestfit_par(mf1$ttp, "aic")$fit$res[,1], mf1$ttp[[best_ttp]]$result$res[,1])
  expect_equal(find_bestfit_par(mf1$pfs, "aic")$fit$res[,1], mf1$pfs[[best_pfs]]$result$res[,1])
  expect_equal(find_bestfit_par(mf1$os, "aic")$fit$res[,1], mf1$os[[best_os]]$result$res[,1])
})

# BIC

ll_ttp <- 1:Ndists |> purrr::map_vec(~mf1$ttp[[.x]]$result$loglik)
ll_pfs <- 1:Ndists |> purrr::map_vec(~mf1$pfs[[.x]]$result$loglik)
ll_os <- 1:Ndists |> purrr::map_vec(~mf1$os[[.x]]$result$loglik)
np_ttp <- 1:Ndists |> purrr::map_vec(~mf1$ttp[[.x]]$result$npars)
np_pfs <- 1:Ndists |> purrr::map_vec(~mf1$pfs[[.x]]$result$npars)
np_os <- 1:Ndists |> purrr::map_vec(~mf1$os[[.x]]$result$npars)
N_ttp <- 1:Ndists |> purrr::map_vec(~mf1$ttp[[.x]]$result$N)
N_pfs <- 1:Ndists |> purrr::map_vec(~mf1$pfs[[.x]]$result$N)
N_os <- 1:Ndists |> purrr::map_vec(~mf1$os[[.x]]$result$N)
bic_ttp <- np_ttp*log(N_ttp)-2*ll_ttp
bic_pfs <- np_pfs*log(N_pfs)-2*ll_pfs
bic_os <- np_os*log(N_os)-2*ll_os
bestb_ttp <- which.min(bic_ttp)
bestb_pfs <- which.min(bic_pfs)
bestb_os <- which.min(bic_os)

test_that("findbest_survreg finds the best fits by min BIC", {
  expect_equal(min(find_bestfit_par(mf1$ttp, "bic")$results$bic), bic_ttp[bestb_ttp])
  expect_equal(min(find_bestfit_par(mf1$pfs, "bic")$results$bic), bic_pfs[bestb_pfs])
  expect_equal(min(find_bestfit_par(mf1$os, "bic")$results$bic), bic_os[bestb_os])
  expect_equal(find_bestfit_par(mf1$ttp, "bic")$fit$res[,1], mf1$ttp[[bestb_ttp]]$result$res[,1])
  expect_equal(find_bestfit_par(mf1$pfs, "bic")$fit$res[,1], mf1$pfs[[bestb_pfs]]$result$res[,1])
  expect_equal(find_bestfit_par(mf1$os, "bic")$fit$res[,1], mf1$os[[bestb_os]]$result$res[,1])
})

