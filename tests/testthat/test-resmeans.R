# Testing for resmeans.R
# ====================

bosonc <- create_dummydata("flexbosms")
days_in_year <- 365.25
days_in_week <- 7
weeks_in_year <- days_in_year/days_in_week

alldists <- c("exp", "weibullPH", "llogis", "lnorm", "gamma",
              "gompertz", "gengamma")


# Fit all distributions to all endpoints
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

# Call the RMD functions
rmd_all <- calc_allrmds(bosonc,
            cuttime = 0,
            Ty = 10,
            dpam = params)

# Check that intermediate functions work
# --------------------------------------

# Calculation of RMDs by model and state
#   rmd_pf_stm: PF state, according to STM-CR or STM-CF models
#   rmd_pd_stm_cr: PD state, according to STM-CR model
#   rmd_pd_stm_cf: PD state, according to STM-CF model
#   rmd_pf_psm: PF state, according to PSM model
#   rmd_os_psm: OS states, according to PSM model

# Check PSM calculations

Ty <- 10
Tw <- Ty * weeks_in_year
# params$pfs$dlist$name is "exp"
exp_psm_pf <- flexsurv::rmst_exp(Tw, rate=params$pfs$res[1])
# params$os$dlist$name is "weibullPH"
exp_psm_os <- flexsurv::rmst_weibullPH(Tw, shape=params$os$res[1,1], scale=params$os$res[2,1])
exp_psm_pd <- exp_psm_os - exp_psm_pf

test_that("PSM results match expected durations", {
  expect_equal(as.numeric(rmd_all$results$pf[1]),
               as.numeric(exp_psm_pf)
               )
  expect_equal(as.numeric(rmd_pf_psm(params, Ty=Ty)),
               as.numeric(exp_psm_pf)
               )
  expect_equal(as.numeric(rmd_all$results$pd[1]),
               as.numeric(exp_psm_pd)
              )
  # rmd_pd_psm does not exist
  expect_equal(as.numeric(rmd_all$results$os[1]),
               as.numeric(exp_psm_os)
               )
  expect_equal(as.numeric(rmd_os_psm(params, Ty=Ty)),
               as.numeric(exp_psm_os)
              )
})

# Check STM-CF calculations

ppd.spec <- list(dist=params$ppd$dlist$name, pars=params$ppd$res[,1])
ttp.spec <- list(dist=params$ttp$dlist$name, pars=params$ttp$res[,1])
pps_cf.spec <- list(dist=params$pps_cf$dlist$name, pars=params$pps_cf$res[,1])

int1 <- function(x) {
  sttp <- calc_surv(x, "par", ttp.spec)
  sppd <- calc_surv(x, "par", ppd.spec)
  sttp*sppd
}
exp_stmcf_pf <- stats::integrate(int1, 0, Tw)$value
exp_stmcr_pf <- exp_stmcf_pf

int2 <- function(x) {
  sttp <- calc_surv(x[1], "par", ttp.spec)
  sppd <- calc_surv(x[1], "par", ppd.spec)
  http <- calc_haz(x[1], "par", ttp.spec)
  sos1 <- calc_surv(x[1], "par", pps_cf.spec)
  sos2 <- calc_surv(x[2], "par", pps_cf.spec)
  ifelse(sos1==0, 0, sttp*sppd*http*sos2/sos1)
}
S <- cbind(c(0,0),c(0, Tw),c(Tw, Tw))
exp_stmcf_pd <- SimplicialCubature::adaptIntegrateSimplex(int2, S)$integral
exp_stmcf_os <- exp_stmcf_pf + exp_stmcf_pd

test_that("STM-CF results match expected durations", {
  expect_equal(as.numeric(rmd_all$results$pf[2]),
               as.numeric(exp_stmcf_pf)
  )
  expect_equal(as.numeric(rmd_pf_stm(params, Ty=Ty)),
               as.numeric(exp_stmcf_pf)
  )
  expect_equal(as.numeric(rmd_all$results$pd[2]),
               as.numeric(exp_stmcf_pd)
  )
  expect_equal(as.numeric(rmd_pd_stm_cf(params, Ty=Ty)),
               as.numeric(exp_stmcf_pd)
  )
  expect_equal(as.numeric(rmd_all$results$os[2]),
               as.numeric(exp_stmcf_os)
  )
  # No rmd_os_stm_cf function
})

# Check STM-CR calculations

pps_cr.spec <- list(dist=params$pps_cr$dlist$name, pars=params$pps_cr$res[,1])

int3 <- function(x) {
  sttp <- calc_surv(x[1], "par", ttp.spec)
  sppd <- calc_surv(x[1], "par", ppd.spec)
  http <- calc_haz(x[1], "par", ttp.spec)
  spps <- calc_surv(x[2]-x[1], "par", pps_cr.spec)
  sttp*sppd*http*spps
}

exp_stmcr_pd <- SimplicialCubature::adaptIntegrateSimplex(int3, S)$integral
exp_stmcr_os <- exp_stmcr_pf + exp_stmcr_pd

test_that("STM-CR results match expected durations", {
  expect_equal(as.numeric(rmd_all$results$pf[3]),
               as.numeric(exp_stmcr_pf)
  )
  expect_equal(as.numeric(rmd_pf_stm(params, Ty=Ty)),
               as.numeric(exp_stmcr_pf)
  )
  expect_equal(as.numeric(rmd_all$results$pd[3]),
               as.numeric(exp_stmcr_pd)
  )
  expect_equal(as.numeric(rmd_pd_stm_cr(params, Ty=Ty)),
               as.numeric(exp_stmcr_pd)
  )
  expect_equal(as.numeric(rmd_all$results$os[3]),
               as.numeric(exp_stmcr_os)
  )
 # No rmd_os_stm_cr function
})

# Expected results
Ty <- 10
Tw <- Ty * weeks_in_year
exp_stmcf_pf2 <- stats::integrate(int1, 0, 10*weeks_in_year)$value
S <- cbind(c(0,0),c(0, Tw),c(Tw, Tw))
exp_stmcr_pd2 <- SimplicialCubature::adaptIntegrateSimplex(int3, S)$integral
exp_stmcf_pd2 <- SimplicialCubature::adaptIntegrateSimplex(int2, S)$integral
exp_psm_pf2 <- flexsurv::rmst_exp(Tw, rate=params$pfs$res[1]) # Exp
exp_psm_os2 <- flexsurv::rmst_weibullPH(Tw,
                             shape=params$os$res[1,1],
                             scale=params$os$res[2,1]) # WeibullPH

test_that("Intermediate results match expected durations", {
  expect_equal(rmd_pf_stm(params, Ty=10), exp_stmcf_pf2)
  expect_equal(rmd_pd_stm_cr(params, Ty=10), exp_stmcr_pd2)
  expect_equal(rmd_pd_stm_cf(params, Ty=10), exp_stmcf_pd2)
  expect_equal(rmd_pf_psm(params, Ty=10), exp_psm_pf2)
  expect_equal(rmd_os_psm(params, Ty=10), exp_psm_os2)
})

test_that("Safe functions produce the same values as originals", {
  expect_equal(rmd_pf_stm(params, Ty=15), prmd_pf_stm(params, Ty=15))
  expect_equal(rmd_pd_stm_cr(params, Ty=15), prmd_pd_stm_cr(params, Ty=15))
  expect_equal(rmd_pd_stm_cf(params, Ty=15), prmd_pd_stm_cf(params, Ty=15))
  expect_equal(rmd_pf_psm(params, Ty=15), prmd_pf_psm(params, Ty=15))
  expect_equal(rmd_os_psm(params, Ty=15), prmd_os_psm(params, Ty=15))
})

Ty_lo <- 0.7
Ty_hi <- 1.1
test_that("Shorter time horizon gives lower mean", {
  expect_lte(rmd_pf_stm(params, Ty=Ty_lo), rmd_pf_stm(params, Ty=Ty_hi))
  expect_lte(rmd_pd_stm_cr(params, Ty=Ty_lo), rmd_pd_stm_cr(params, Ty=Ty_hi))
  expect_lte(rmd_pd_stm_cf(params, Ty=Ty_lo), rmd_pd_stm_cf(params, Ty=Ty_hi))
  expect_lte(rmd_pf_psm(params, Ty=Ty_lo), rmd_pf_psm(params, Ty=Ty_hi))
  expect_lte(rmd_os_psm(params, Ty=Ty_lo), rmd_os_psm(params, Ty=Ty_hi))
})

test_that("Zero time horizon gives zero mean", {
  expect_equal(rmd_pf_stm(params, Ty=0), 0)
  expect_equal(rmd_pd_stm_cr(params, Ty=0), 0)
  expect_equal(rmd_pd_stm_cf(params, Ty=0), 0)
  expect_equal(rmd_pf_psm(params, Ty=0), 0)
  expect_equal(rmd_os_psm(params, Ty=0), 0)
})

# Check discounting reduces the means

test_that("Including discounting reduces the mean", {
  expect_gte(rmd_pf_stm(params, Ty=15), rmd_pf_stm(params, 15, discrate=0.03))
  expect_gte(rmd_pd_stm_cr(params, Ty=15), rmd_pd_stm_cr(params, Ty=15, discrate=0.03))
  expect_gte(rmd_pd_stm_cf(params, Ty=15), rmd_pd_stm_cf(params, Ty=15, discrate=0.03))
  expect_gte(rmd_pf_psm(params, Ty=15), rmd_pf_psm(params, Ty=15, discrate=0.03))
  expect_gte(rmd_os_psm(params, Ty=15), rmd_os_psm(params, Ty=15, discrate=0.03))
})