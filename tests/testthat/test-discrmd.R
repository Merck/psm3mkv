# Tests for discrmd.R
# -------------------

# Create dataset
bosonc <- create_dummydata("flexbosms")

# Multiply durations by 300 so they can be constrained by lifetables
bosonc$pfs.durn <- bosonc$pfs.durn * 300
bosonc$os.durn <- bosonc$os.durn * 300
bosonc$ttp.durn <- bosonc$ttp.durn * 300

# Fit distributions
fits <- fit_ends_mods_spl(bosonc)
params <- list(
 ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
 ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
 pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
 os = find_bestfit_spl(fits$os, "aic")$fit,
 pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
 pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
)

# Add a lifetable constraint
ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)

# Integral results
pf_psm <- prmd_pf_psm(params)
os_psm <- prmd_os_psm(params)
pf_stm <- prmd_pf_stm(params)
pd_stmcf <- prmd_pd_stm_cf(params)
pd_stmcr <- prmd_pd_stm_cr(params)

# Discretized results without lifetables
psm_drmd_wo <- drmd_psm(ptdata=bosonc, dpam=params)
stmcf_drmd_wo <- drmd_stm_cf(dpam=params)
stmcr_drmd_wo <- drmd_stm_cr(dpam=params)

# Discretized results with lifetables
psm_drmd_wi <- drmd_psm(ptdata=bosonc, dpam=params, lifetable=ltable)
stmcf_drmd_wi <- drmd_stm_cf(dpam=params, lifetable=ltable)
stmcr_drmd_wi <- drmd_stm_cr(dpam=params, lifetable=ltable)

# Check that discretized results without lifetables are close to integral results
# 'Close to' = within +/-5%
margin <- 0.05

test_that("Discretized results <= integral results + margin", {
  expect_lte(as.numeric(psm_drmd_wo$pf),
               as.numeric(pf_psm)*(1+margin)
  )
  expect_lte(as.numeric(psm_drmd_wo$os),
             as.numeric(os_psm)*(1+margin)
  )
  expect_lte(as.numeric(psm_drmd_wo$pd),
             (as.numeric(os_psm)-as.numeric(pf_psm))*(1+margin)
  )
  expect_lte(as.numeric(stmcf_drmd_wo$pf),
             as.numeric(pf_stm)*(1+margin)
  )
  expect_lte(as.numeric(stmcf_drmd_wo$os),
             (as.numeric(pf_stm)+as.numeric(pd_stmcf))*(1+margin)
  )
  expect_lte(as.numeric(stmcf_drmd_wo$pd),
             as.numeric(pd_stmcf)*(1+margin)
  )
  expect_lte(as.numeric(stmcr_drmd_wo$pf),
             as.numeric(pf_stm)*(1+margin)
  )
  expect_lte(as.numeric(stmcr_drmd_wo$os),
             (as.numeric(pf_stm)+as.numeric(pd_stmcr))*(1+margin)
  )
  expect_lte(as.numeric(stmcr_drmd_wo$pd),
             as.numeric(pd_stmcr)*(1+margin)
  )
})

test_that("Discretized results >= integral results - margin", {
  expect_gte(as.numeric(psm_drmd_wo$pf),
             as.numeric(pf_psm)/(1+margin)
  )
  expect_gte(as.numeric(psm_drmd_wo$os),
             as.numeric(os_psm)/(1+margin)
  )
  expect_gte(as.numeric(psm_drmd_wo$pd),
             (as.numeric(os_psm)-as.numeric(pf_psm))/(1+margin)
  )
  expect_gte(as.numeric(stmcf_drmd_wo$pf),
             as.numeric(pf_stm)/(1+margin)
  )
  expect_gte(as.numeric(stmcf_drmd_wo$os),
             (as.numeric(pf_stm)+as.numeric(pd_stmcf))/(1+margin)
  )
  expect_gte(as.numeric(stmcf_drmd_wo$pd),
             as.numeric(pd_stmcf)/(1+margin)
  )
  expect_gte(as.numeric(stmcr_drmd_wo$pf),
             as.numeric(pf_stm)/(1+margin)
  )
  expect_gte(as.numeric(stmcr_drmd_wo$os),
             (as.numeric(pf_stm)+as.numeric(pd_stmcr))/(1+margin)
  )
  expect_gte(as.numeric(stmcr_drmd_wo$pd),
             as.numeric(pd_stmcr)/(1+margin)
  )
})

# Check that constraining by a lifetable reduces RMD values

test_that("Discretized results with constraint <= without", {
  expect_lte(as.numeric(psm_drmd_wi$pf),
             as.numeric(psm_drmd_wo$pf)
  )
  expect_lte(as.numeric(psm_drmd_wi$pd),
             as.numeric(psm_drmd_wo$pd)
  )
  expect_lte(as.numeric(psm_drmd_wi$os),
             as.numeric(psm_drmd_wo$os)
  )
  expect_lte(as.numeric(stmcf_drmd_wi$pf),
             as.numeric(stmcf_drmd_wo$pf)
  )
  expect_lte(as.numeric(stmcf_drmd_wi$pd),
             as.numeric(stmcf_drmd_wo$pd)
  )
  expect_lte(as.numeric(stmcf_drmd_wi$os),
             as.numeric(stmcf_drmd_wo$os)
  )
  expect_lte(as.numeric(stmcr_drmd_wi$pf),
             as.numeric(stmcr_drmd_wo$pf)
  )
  expect_lte(as.numeric(stmcr_drmd_wi$pd),
             as.numeric(stmcr_drmd_wo$pd)
  )
  expect_lte(as.numeric(stmcr_drmd_wi$os),
             as.numeric(stmcr_drmd_wo$os)
  )
})