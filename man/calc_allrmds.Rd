% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resmeans.R
\name{calc_allrmds}
\alias{calc_allrmds}
\title{Calculate restricted mean durations for each health state and all three models}
\usage{
calc_allrmds(
  simdat,
  inclset = 0,
  dpam,
  cuttime = 0,
  Ty = 10,
  lifetable = NA,
  discrate = 0
)
}
\arguments{
\item{simdat}{Dataset of patient level data. Must be a tibble with columns named:
\itemize{
\item ptid: patient identifier
\item pfs.durn: duration of PFS from baseline
\item pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
\item os.durn: duration of OS from baseline
\item os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
\item ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
\item ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
}}

\item{inclset}{Vector to indicate which patients to include in analysis}

\item{dpam}{List of statistical fits to each endpoint required in PSM, STM-CF and STM-CR models.}

\item{cuttime}{Time cutoff - this is nonzero for two-piece models.}

\item{Ty}{Time duration over which to calculate. Assumes input is in years, and patient-level data is recorded in weeks.}

\item{lifetable}{Optional, a life table. Columns must include \code{lttime} (time in years, or 52.18 times shorter than the time index elsewhere, starting from zero) and \code{lx}}

\item{discrate}{Discount rate (\% per year)}
}
\value{
List of detailed numeric results
\itemize{
\item cutadj indicates the survival function and area under the curves for PFS and OS up to the cutpoint
\item results provides results of the restricted means calculations, by model and state.
}
}
\description{
Calculate restricted mean durations for each health state (progression free and progressed disease) for all three models (partitioned survival, clock forward state transition model, clock reset state transition model).
}
\examples{
# Create dataset and fit survival models (splines)
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
# Find mean(s)
calc_allrmds(bosonc, dpam=params)
}
\seealso{
Restricted means are provided by \code{\link[=rmd_pf_psm]{rmd_pf_psm()}}, \code{\link[=rmd_os_psm]{rmd_os_psm()}}, \code{\link[=rmd_pf_stm]{rmd_pf_stm()}}, \code{\link[=rmd_pd_stm_cf]{rmd_pd_stm_cf()}} and \code{\link[=rmd_pd_stm_cr]{rmd_pd_stm_cr()}}. The function \link{calc_allrmds_boot} provides a version for bootstrapping.
}
