# Testing of functions in basics.R
# Functions are: calc_pdist, calc_surv, calc_haz, gen_rand, calc_rmd, give_noparams
# ==============================================================================

# 1. give_noparams
# ----------------

test_that("each distribution works", {
  expect_equal(psm3mkv::give_noparams_par("exp"), 1)
  expect_equal(psm3mkv::give_noparams_par("weibullPH"), 2)
  expect_equal(psm3mkv::give_noparams_par("weibull"), 2)
  expect_equal(psm3mkv::give_noparams_par("llogis"), 2)
  expect_equal(psm3mkv::give_noparams_par("lnorm"), 2)
  expect_equal(psm3mkv::give_noparams_par("gamma"), 2)
  expect_equal(psm3mkv::give_noparams_par("gompertz"), 2)
  expect_equal(psm3mkv::give_noparams_par("gengamma"), 3)
  expect_equal(psm3mkv::give_noparams_par("gengamma.orig"), 3)
})

test_that("entering >1 distribution yields an error", {
  expect_error(psm3mkv::give_noparams_par(c("exp", "lnorm")))
  expect_error(psm3mkv::give_noparams_par(c("gengamma", "turnip")))
  expect_error(psm3mkv::give_noparams_par(c("carrot", "llogis")))
})

# Check with splines

# Three example spline fits
fit_spl1 <- flexsurv::flexsurvspline(
  survival::Surv(recyrs, censrec) ~ 1,
  data=flexsurv::bc,
  k=1,
  scale="odds")
fit_spl2 <- flexsurv::flexsurvspline(
  survival::Surv(recyrs, censrec) ~ 1,
  data=flexsurv::bc,
  k=2,
  scale="hazard")
fit_spl3 <- flexsurv::flexsurvspline(
  survival::Surv(recyrs, censrec) ~ 1,
  data=flexsurv::bc,
  k=3,
  scale="normal")

# Three example list specifications
spec_spl1 <- list(gammas = fit_spl1$coefficients,
             knots = fit_spl1$aux$knots,
             scale = fit_spl1$scale)
spec_spl2 <- list(gammas = fit_spl2$coefficients,
                 knots = fit_spl2$aux$knots,
                 scale = fit_spl2$scale)
spec_spl3 <- list(gammas = fit_spl3$coefficients,
                 knots = fit_spl3$aux$knots,
                 scale = fit_spl3$scale)

test_that("Spline specifications give correct parameter numbers", {
  expect_equal(psm3mkv::give_noparams(type="Splines", spec=spec_spl1), 3)
  expect_equal(psm3mkv::give_noparams(type="Splines", spec=spec_spl2), 4)
  expect_equal(psm3mkv::give_noparams(type="Splines", spec=spec_spl3), 5)
})

test_that("Parametric specification give correct parameter numbers", {
  expect_equal(psm3mkv::give_noparams(type="par", spec=list(dist="exp")), 1)
  expect_equal(psm3mkv::give_noparams(type="Parametrci", spec=list(dist="weibullPH")), 2)
  expect_equal(psm3mkv::give_noparams(type="paR", spec=list(dist="gengamma")), 3)
})

# 2. calc_rmd
# -----------

test_that("Restricted mean equals mean over inf horizon, param", {
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="exp", pars=c(0.2))),
    flexsurv::mean_exp(0.2)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="weibullPH", pars=c(1, 1))),
    flexsurv::mean_weibullPH(1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="weibull", pars=c(1, 1))),
    flexsurv::mean_weibull(1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="llogis", pars=c(2, 3))),
    flexsurv::mean_llogis(2, 3)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="lnorm", pars=c(1,1))),
    flexsurv::mean_lnorm(1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="gamma", pars=c(1, 1))),
    flexsurv::mean_gamma(1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="gompertz", pars=c(1, 1))),
    flexsurv::mean_gompertz(1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="gengamma", pars=c(1, 1, 1))),
    flexsurv::mean_gengamma(1, 1, 1)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="Parametric", spec=list(dist="gengamma.orig", pars=c(1, 1, 1))),
    flexsurv::mean_gengamma.orig(1, 1, 1)
    )
})

test_that("Restricted mean equals mean over inf horizon, splines", {
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="splines", spec=spec_spl1),
    flexsurv::mean_survspline(gamma=spec_spl1$gammas,
                              knots=spec_spl1$knots,
                              scale=spec_spl1$scale)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="splines", spec=spec_spl2),
    flexsurv::mean_survspline(gamma=spec_spl2$gammas,
                              knots=spec_spl2$knots,
                              scale=spec_spl2$scale)
    )
  expect_equal(
    psm3mkv::calc_rmd(Tw=Inf, type="splines", spec=spec_spl3),
    flexsurv::mean_survspline(gamma=spec_spl3$gammas,
                              knots=spec_spl3$knots,
                              scale=spec_spl3$scale)
    )
})

test_that("Restricted mean < mean, parametric", {
  expect_lt(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="exp", pars=c(0.2))),
    flexsurv::mean_exp(0.2)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=5, type="Parametric", spec=list(dist="weibullPH", pars=c(1, 1))),
    flexsurv::mean_weibullPH(1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=15, type="Parametric", spec=list(dist="weibull", pars=c(1, 1))),
    flexsurv::mean_weibull(1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=50, type="Parametric", spec=list(dist="llogis", pars=c(2, 3))),
    flexsurv::mean_llogis(2, 3)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=200, type="Parametric", spec=list(dist="lnorm", pars=c(1,1))),
    flexsurv::mean_lnorm(1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=0.1, type="Parametric", spec=list(dist="gamma", pars=c(1, 1))),
    flexsurv::mean_gamma(1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=1, type="Parametric", spec=list(dist="gompertz", pars=c(1, 1))),
    flexsurv::mean_gompertz(1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=20, type="Parametric", spec=list(dist="gengamma", pars=c(1, 1, 1))),
    flexsurv::mean_gengamma(1, 1, 1)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="gengamma.orig", pars=c(1, 1, 1))),
    flexsurv::mean_gengamma.orig(1, 1, 1)
  )
})

test_that("Restricted mean < mean, splines", {
  expect_lt(
    psm3mkv::calc_rmd(Tw=10, type="splines", spec=spec_spl1),
    flexsurv::mean_survspline(gamma=spec_spl1$gammas,
                              knots=spec_spl1$knots,
                              scale=spec_spl1$scale)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=15, type="splines", spec=spec_spl2),
    flexsurv::mean_survspline(gamma=spec_spl2$gammas,
                              knots=spec_spl2$knots,
                              scale=spec_spl2$scale)
  )
  expect_lt(
    psm3mkv::calc_rmd(Tw=30, type="splines", spec=spec_spl3),
    flexsurv::mean_survspline(gamma=spec_spl3$gammas,
                              knots=spec_spl3$knots,
                              scale=spec_spl3$scale)
  )
})

test_that("Calling restricted mean function correctly, parametric", {
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="exp", pars=0.2)),
    flexsurv::rmst_exp(t=10, rate=0.2, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="weibullPH", pars=c(1, 1))),
    flexsurv::rmst_weibullPH(t=10, shape=1, scale=1, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="weibull", pars=c(1, 1))),
    flexsurv::rmst_weibull(t=10, shape=1, scale=1, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="weibullPH", pars=c(1, 1))),
    flexsurv::rmst_weibullPH(t=10, shape=1, scale=1, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="llogis", pars=c(1, 1))),
    flexsurv::rmst_llogis(t=10, shape=1, scale=1, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=10, type="Parametric", spec=list(dist="lnorm", pars=c(4, 3))),
    flexsurv::rmst_lnorm(t=10, meanlog=4, sdlog=3, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=20, type="Parametric", spec=list(dist="gamma", pars=c(3, 2))),
    flexsurv::rmst_gamma(t=20, shape=3, rate=2, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=20, type="Parametric", spec=list(dist="gompertz", pars=c(3, 2))),
    flexsurv::rmst_gompertz(t=20, shape=3, rate=2, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=15, type="Parametric", spec=list(dist="gengamma", pars=c(3, 2, 1))),
    flexsurv::rmst_gengamma(t=15, mu=3, sigma=2, Q=1, start=0)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=15, type="Parametric", spec=list(dist="gengamma.orig", pars=c(2, 3, 1))),
    flexsurv::rmst_gengamma.orig(t=15, shape=2, scale=3, k=1, start=0)
  )
})

test_that("Calling restricted mean function correctly, splines", {
  expect_equal(
    psm3mkv::calc_rmd(Tw=30, type="splines", spec=spec_spl1),
    flexsurv::rmst_survspline(t=30,
                          gamma=spec_spl1$gammas,
                          knots=spec_spl1$knots,
                          scale=spec_spl1$scale)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=15, type="splines", spec=spec_spl2),
    flexsurv::rmst_survspline(t=15,
                          gamma=spec_spl2$gammas,
                          knots=spec_spl2$knots,
                          scale=spec_spl2$scale)
  )
  expect_equal(
    psm3mkv::calc_rmd(Tw=20, type="splines", spec=spec_spl3),
    flexsurv::rmst_survspline(t=20,
                          gamma=spec_spl3$gammas,
                          knots=spec_spl3$knots,
                          scale=spec_spl3$scale)
  )
})

# 3. calc_surv and gen_rand
# -------------------------

# calc_surv(time=0) = 0

test_that("survival at time zero is 1, params", {
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="exp", pars=0.2)),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="weibullPH", pars=c(1,1))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="weibull", pars=c(1,1))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="llogis", pars=c(4,3))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="lnorm", pars=c(2,3))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="gamma", pars=c(2,1))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="gompertz", pars=c(0.3,0.01))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="gengamma", pars=c(2.5,1.5,0.5))),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="par", spec=list(dist="gengamma.orig", pars=c(0.1,10,0.5))),
    1
  )
})

test_that("survival at time zero is 1, splines", {
  expect_equal(
    psm3mkv::calc_surv(time=0, type="splines", spec=spec_spl1),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="splines", spec=spec_spl2),
    1
  )
  expect_equal(
    psm3mkv::calc_surv(time=0, type="splines", spec=spec_spl3),
    1
  )
})

