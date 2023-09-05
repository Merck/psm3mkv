% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probgraphs.R
\name{prob_pf_stm}
\alias{prob_pf_stm}
\title{Calculate probability of being progression free in either markov model (clock forward or clock reset)}
\usage{
prob_pf_stm(time, dpam, starting = c(1, 0, 0))
}
\arguments{
\item{time}{Time (numeric and vectorized)}

\item{dpam}{List of survival regressions for model endpoints. This must include pre-progression death (PPD) and time to progression (TTP).}

\item{starting}{Vector of membership probabilities (PF, PD, death) at time zero.}
}
\value{
Numeric value
}
\description{
Calculates membership probability for the progression free state, at a particular time (vectorized), given either markov model (clock forward or clock reset) with given statistical distributions and parameters.
}
\examples{
bosonc <- create_dummydata("flexbosms")
fits <- fit_ends_mods_spl(bosonc)
# Pick out best distribution according to min AIC
params <- list(
  ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
  ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
  pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
  os = find_bestfit_spl(fits$os, "aic")$fit,
  pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
  pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
)
prob_pf_stm(0:100, params)
}