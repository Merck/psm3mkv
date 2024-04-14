# Tests for ltablesurv.R, which contains
# - vlookup()
# - calc_ltsurv()
# - calc_ex() - which is not actually used anywhere (yet) so not tested

# Vlookup is designed to work like Excel's vlookup() function, except with a choice of interpolation methods

# First let's test single value lookups without any interpolation

# Example lookup-table
looktab1 <- tibble::tibble(
  rowno=1:10,
  ind=(rowno-6)*10,
  val=runif(10)
)

test_that("Exact lookups work", {
  expect_equal(vlookup(looktab1$ind[3], looktab1$ind, looktab1$val), looktab1$val[3])
  expect_equal(vlookup(looktab1$ind[1], looktab1$ind, looktab1$val), looktab1$val[1])
  expect_equal(vlookup(looktab1$ind[10], looktab1$ind, looktab1$val), looktab1$val[10])
})

test_that("Error occurs if lookup value out of range", {
  expect_error(vlookup(-90, looktab1$ind, looktab1$val))
  expect_error(vlookup(190, looktab1$ind, looktab1$val))
})

test_that("All methods produce the same result, with exact lookup", {
  expect_equal(vlookup(looktab1$ind[3], looktab1$ind, looktab1$val, method="floor"), looktab1$val[3])
  expect_equal(vlookup(looktab1$ind[3], looktab1$ind, looktab1$val, method="ceiling"), looktab1$val[3])
  expect_equal(vlookup(looktab1$ind[3], looktab1$ind, looktab1$val, method="arith"), looktab1$val[3])
  expect_equal(vlookup(looktab1$ind[3], looktab1$ind, looktab1$val, method=""), looktab1$val[3])
  expect_equal(vlookup(looktab1$ind[1], looktab1$ind, looktab1$val, method="floor"), looktab1$val[1])
  expect_equal(vlookup(looktab1$ind[1], looktab1$ind, looktab1$val, method="ceiling"), looktab1$val[1])
  expect_equal(vlookup(looktab1$ind[1], looktab1$ind, looktab1$val, method="arith"), looktab1$val[1])
  expect_equal(vlookup(looktab1$ind[1], looktab1$ind, looktab1$val, method="garbage"), looktab1$val[1])
  expect_equal(vlookup(looktab1$ind[10], looktab1$ind, looktab1$val, method="floor"), looktab1$val[10])
  expect_equal(vlookup(looktab1$ind[10], looktab1$ind, looktab1$val, method="ceiling"), looktab1$val[10])
  expect_equal(vlookup(looktab1$ind[10], looktab1$ind, looktab1$val, method="arith"), looktab1$val[10])
  expect_equal(vlookup(looktab1$ind[10], looktab1$ind, looktab1$val, method=""), looktab1$val[10])
})

# Now try multi-value but exact lookups

mval <- c(1, 3, 10)

test_that("Multivalue exact lookups produce the same combined result", {
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method="floor"), looktab1$val[mval])
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method="ceiling"), looktab1$val[mval])
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method="arith"), looktab1$val[mval])
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method="geom"), looktab1$val[mval])
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method="garbage"), looktab1$val[mval])
  expect_equal(vlookup(looktab1$ind[mval], looktab1$ind, looktab1$val, method=NA), looktab1$val[mval])
})

# Single values with interpolations

looktab1$ind[8] <- 18 # was 20
indval1.act <- 24
indval1.lo <- 18
indval1.hi <- 30
lookval1.lo <- looktab1$val[8]
lookval1.hi <- looktab1$val[9]
lookval1.arith <- (lookval1.lo*(indval1.hi-indval1.act) + lookval1.hi*(indval1.act-indval1.lo))/(indval1.hi-indval1.lo)
lookval1.geom <- (lookval1.lo^(indval1.hi-indval1.act) * lookval1.hi^(indval1.act-indval1.lo))^(1/(indval1.hi-indval1.lo))

test_that("Single value interpolation produces the right result with each method", {
  expect_equal(vlookup(indval1.act, looktab1$ind, looktab1$val, method="floor"), lookval1.lo)
  expect_equal(vlookup(indval1.act, looktab1$ind, looktab1$val, method="ceiling"), lookval1.hi)
  expect_equal(vlookup(indval1.act, looktab1$ind, looktab1$val, method="arith"), lookval1.arith)
  expect_equal(vlookup(indval1.act, looktab1$ind, looktab1$val, method="geom"), lookval1.geom)
  expect_equal(vlookup(indval1.act, looktab1$ind, looktab1$val, method=""), lookval1.geom)
})

# Multi values with interpolations

indval2.act <- -36
indval2.lo <- -40
indval2.hi <- -30
lookval2.lo <- looktab1$val[2]
lookval2.hi <- looktab1$val[3]
lookval2.arith <- (lookval2.lo*(indval2.hi-indval2.act) + lookval2.hi*(indval2.act-indval2.lo))/(indval2.hi-indval2.lo)
lookval2.geom <- (lookval2.lo^(indval2.hi-indval2.act) * lookval2.hi^(indval2.act-indval2.lo))^(1/(indval2.hi-indval2.lo))
combind.act <- c(indval1.act, indval2.act)
combval.lo <- c(lookval1.lo, lookval2.lo)
combval.hi <- c(lookval1.hi, lookval2.hi)
combval.arith <- c(lookval1.arith, lookval2.arith)
combval.geom <- c(lookval1.geom, lookval2.geom)

test_that("Multi value interpolation produces the right result with each method", {
  expect_equal(vlookup(combind.act, looktab1$ind, looktab1$val, method="floor"), combval.lo)
  expect_equal(vlookup(combind.act, looktab1$ind, looktab1$val, method="ceiling"), combval.hi)
  expect_equal(vlookup(combind.act, looktab1$ind, looktab1$val, method="arith"), combval.arith)
  expect_equal(vlookup(combind.act, looktab1$ind, looktab1$val, method="geom"), combval.geom)
  expect_equal(vlookup(combind.act, looktab1$ind, looktab1$val, method=""), combval.geom)
})

# Now test the lifetime survival function, which requires a lifetable
# Lifetable must have lttime and lx columns
ltable <- data.frame(lttime=0:20, lx=1000-10*(0:20))

timecheck <- 65
rowcheck <- 6

test_that("Single value survivals are correct", {
  expect_equal(calc_ltsurv(5, ltable), ltable$lx[6]/ltable$lx[1])
  expect_equal(calc_ltsurv(0, ltable), 1)
  expect_equal(calc_ltsurv(20, ltable), ltable$lx[21]/ltable$lx[1])
})

test_that("Return error when lookup is out of range", {
  expect_error(calc_ltsurv(20.1, ltable))
  expect_error(calc_ltsurv(-0.1, ltable))
})

# Now check that the survival functions work with vectors

test_that("Multi value survivals are correct when in range, regardless of order", {
  expect_equal(calc_ltsurv(c(1,2,3), ltable), c(calc_ltsurv(1, ltable), calc_ltsurv(2, ltable), calc_ltsurv(3, ltable)))
  expect_equal(calc_ltsurv(c(4,13,8), ltable), c(calc_ltsurv(4, ltable), calc_ltsurv(13, ltable), calc_ltsurv(8, ltable)))
})

test_that("Multi value survivals return error when any index is out of range", {
  expect_error(calc_ltsurv(c(-5,5,55), ltable))
  expect_error(calc_ltsurv(c(-5,55,5), ltable))
})

# Finally check the density

# Derive density
ltable <- ltable |>
  dplyr::mutate(
    dx = lx - dplyr::lead(lx),
    qx = dx/lx,
    hx = -log(1-qx),
    surv = lx/lx[1],
    dens = hx * surv
  )

test_that("Density function correctly derived with regular lifetable", {
  expect_equal(calc_ltdens(ltable$lttime, ltable), ltable$dens)
  expect_equal(calc_ltdens(0:19+0.2, ltable, method="floor"), ltable$dens[0:20])
})

# Try again with some times missing
ltable <- ltable[c(1,3,5,7,8,10,14,15,17,19),] |>
  dplyr::select(lttime, lx) |>
  dplyr:: mutate(
    dx = lx - dplyr::lead(lx),
    tx = dplyr::lead(lttime) - lttime,
    qx = dx/lx,
    hx = -log(1-qx)/tx,
    surv = lx/lx[1],
    dens = hx * surv
  )

test_that("Density function works with times missing", {
  expect_equal(calc_ltdens(ltable$lttime, ltable), ltable$dens)
  expect_equal(calc_ltdens(ltable$lttime[1:9]+0.2, ltable, method="floor"), ltable$dens[1:9])
})

