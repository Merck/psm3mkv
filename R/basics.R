#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the psm3mkv program.
#
#  psm3mkv is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==================================================================
# Basic helper functions
# All these functions heavily rely on the flexsurv or stats packages
# basics.R
# ==================================================================

#' Calculate the distribution function for parametric functions
#' @description Calculate the value of the distribution function for a parametric distribution, given the statistical distribution and its parameters.
#' @param time is the time at which the distribution function should be calculated.
#' @param dist is the statistical distribution (named per [flexsurv::flexsurvreg]).
#' @param pars is a vector of the parameters for that distribution.
#' - Exponential distribution (`exp`) requires the rate parameter.
#' - Weibull distribution (both `weibullPH` and `weibull` formulations) requires the shape and scale parameters.
#' - Log-logistic distribution (`llogis`) requires the shape and scale parameters.
#' - Log-normal distribution (`lnorm`) requires the meanlog and sdlog parameters.
#' - Gamma and Gompertz distributions (`gamma` and `gompertz`) require the shape and rate parameters.
#' - Generalized Gamma requires the mu, sigma and Q parameters if using the standard parameterization (`gengamma`) or shape, scale and k parameters if using the original parameterization (`gengamma.orig`).
#' @seealso [flexsurv::flexsurvreg]
#' @return The value of the distribution function, a numeric value.
#' @export
#' @examples
#' calc_pdist_par(10, "exp", 0.01)
#' calc_pdist_par(5, "lnorm", c(3, 1))
calc_pdist_par <- function(time, dist, pars) {
  if (dist=="exp") {
    stats::pexp(time, rate=pars[1])
  } else if (dist=="weibullPH") {
    flexsurv::pweibullPH(time, shape=pars[1], scale=pars[2])
  } else if (dist=="weibull" | dist=="weibull.quiet") {
    stats::pweibull(time, shape=pars[1], scale=pars[2])
  } else if (dist=="llogis") {
    flexsurv::pllogis(time, shape=pars[1], scale=pars[2])
  } else if (dist=="lnorm") {
    stats::plnorm(time, meanlog=pars[1], sdlog=pars[2])
  } else if (dist=="gamma") {
    stats::pgamma(time, shape=pars[1], rate=pars[2])
  } else if (dist=="gompertz") {
    flexsurv::pgompertz(time, shape=pars[1], rate=pars[2])
  } else if (dist=="gengamma") {
    flexsurv::pgengamma(time, mu=pars[1], sigma=pars[2], Q=pars[3])
  } else if (dist=="gengamma.orig") {
    flexsurv::pgengamma.orig(time, shape=pars[1], scale=pars[2], k=pars[3])
  } else {
    NA
  }
}

#' Calculate the distribution function
#' @description Calculate the value of the distribution function, given a regular parametric or Royston-Parmar formulation
#' @param time is the time at which the distribution function should be calculated.
#' @param type is either "par" for regular parametric form (exponential, weibull etc) or "spl" for Royston-Parmar splines.
#' @param spec is a list comprising:
#' If type=="par": `dist` is the statistical distribution (named per [flexsurv::flexsurvreg]) and `pars` is a vector of the parameters for that distribution.
#' - Exponential distribution (`exp`) requires the rate parameter.
#' - Weibull distribution (both `weibullPH` and `weibull` formulations) requires the shape and scale parameters.
#' - Log-logistic distribution (`llogis`) requires the shape and scale parameters.
#' - Log-normal distribution (`lnorm`) requires the meanlog and sdlog parameters.
#' - Gamma and Gompertz distributions (`gamma` and `gompertz`) require the shape and rate parameters.
#' - Generalized Gamma requires the mu, sigma and Q parameters if using the standard parameterization (`gengamma`) or shape, scale and k parameters if using the original parameterization (`gengamma.orig`).
#' If type=="spl":
#' - `gamma` - Vector of parameters describing the baseline spline function, as described in [flexsurv::flexsurvspline]. This may be supplied as a vector with number of elements equal to the length of knots, in which case the parameters are common to all times. Alternatively a matrix may be supplied, with rows corresponding to different times, and columns corresponding to knots.
#' - `knots` - Vector of locations of knots on the axis of log time, supplied in increasing order. Unlike in [flexsurv::flexsurvspline], these include the two boundary knots.
#' - `scale` - Either "hazard", "odds", or "normal", as described in [flexsurv::flexsurvspline]. With the default of no knots in addition to the boundaries, this model reduces to the Weibull, log-logistic and log-normal respectively. The scale must be common to all times.
#' @seealso [flexsurv::flexsurvspline] and [flexsurv::flexsurvreg]
#' @return The value of the distribution function, a numeric value.
#' @export
#' @examples
#' calc_pdist(time=1:5,
#'     type="spl",
#'     spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
#'     )
#' calc_pdist(time=1:5,
#'     type="par",
#'     spec=list(dist="lnorm", pars=c(3,1))
#'     )
calc_pdist <- function(time, type, spec){
  type <- tolower(substr(type, 1, 3))
  if (type=="spl") {
    flexsurv::psurvspline(q = time,
                          gamma = spec$gamma,
                          beta = 0,
                          X = 0,
                          knots = spec$knots,
                          scale = spec$scale,
                          timescale = "log",
                          offset = 0,
                          lower.tail = TRUE,
                          log.p = FALSE)
  } else if (type=="par") {
    calc_pdist_par(time, dist=spec$dist, pars=spec$pars)
  } else {
    NA
  }
}

#' Calculate the value of the survival function
#' @description Calculate the value of the survival function, given the statistical distribution and its parameters.
#' @inheritParams calc_pdist
#' @return The value of the survival function, a numeric value.
#' @inherit calc_pdist seealso
#' @export
#' @examples
#' calc_surv(time=1:5,
#'     type="spl",
#'     spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
#'     )
#' calc_surv(time=1:5,
#'     type="par",
#'     spec=list(dist="lnorm", pars=c(3,1))
#'     )
calc_surv <- function(time, type, spec) {
  1-calc_pdist(time, type, spec)
}

#' Calculate the value of the hazard function (parametric form)
#' @description Calculate the value of the hazard function, given the statistical distribution and its parameters.
#' @inheritParams calc_pdist_par
#' @inherit calc_pdist_par seealso
#' @return The value of a hazard function, a numeric value.
#' @export
#' @examples
#' calc_haz_par(10, "exp", 0.01)
#' calc_haz_par(5, "lnorm", c(3, 1))
calc_haz_par <- function(time, dist, pars) {
  if (dist=="exp") {
    flexsurv::hexp(time, rate=pars[1])
  } else if (dist=="weibullPH") {
    flexsurv::hweibullPH(time, shape=pars[1], scale=pars[2])
  } else if (dist=="weibull" | dist=="weibull.quiet") {
    flexsurv::hweibull(time, shape=pars[1], scale=pars[2])
  } else if (dist=="llogis") {
    flexsurv::hllogis(time, shape=pars[1], scale=pars[2])
  } else if (dist=="lnorm") {
    flexsurv::hlnorm(time, meanlog=pars[1], sdlog=pars[2])
  } else if(dist=="gamma") {
    flexsurv::hgamma(time, shape=pars[1], rate=pars[2])
  } else if(dist=="gompertz") {
    flexsurv::hgompertz(time, shape=pars[1], rate=pars[2])
  } else if(dist=="gengamma") {
    flexsurv::hgengamma(time, mu=pars[1], sigma=pars[2], Q=pars[3])
  } else if(dist=="gengamma.orig") {
    flexsurv::hgengamma.orig(time, shape=pars[1], scale=pars[2], k=pars[3])
  } else {
    NA
  }
}

#' Calculate the value of the hazard function
#' @description Calculate the value of the hazard function, given specification as either parametric or Royston-Parmar splines model
#' @inheritParams calc_surv
#' @inherit calc_haz_par return
#' @export
#' @examples
#' calc_haz(time=1:5,
#'     type="spl",
#'     spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
#'     )
#' calc_haz(time=1:5,
#'     type="par",
#'     spec=list(dist="lnorm", pars=c(3,1))
#'     )
calc_haz <- function(time, type, spec){
  type <- tolower(substr(type, 1, 3))
  if (type=="spl") {
    flexsurv::hsurvspline(x = time,
                                 gamma = spec$gamma,
                                 beta = 0,
                                 X = 0,
                                 knots = spec$knots,
                                 scale = spec$scale,
                                 timescale = "log",
                                 offset = 0)
  } else if (type=="par") {
    calc_haz_par(time, dist=spec$dist, pars=spec$pars)
  } else {
    NA
  }
}

#' Calculate the value of the density function (parametric form)
#' @description Calculate the value of the density function, given the statistical distribution and its parameters.
#' @inheritParams calc_pdist_par
#' @inherit calc_pdist_par seealso
#' @return The value of a density function, a numeric value.
#' @export
#' @examples
#' calc_dens_par(10, "exp", 0.01)
#' calc_dens_par(5, "lnorm", c(3, 1))
calc_dens_par <- function(time, dist, pars) {
  if (dist=="exp") {
    stats::dexp(time, rate=pars[1])
  } else if (dist=="weibullPH") {
    flexsurv::dweibullPH(time, shape=pars[1], scale=pars[2])
  } else if (dist=="weibull" | dist=="weibull.quiet") {
    stats::dweibull(time, shape=pars[1], scale=pars[2])
  } else if (dist=="llogis") {
    flexsurv::dllogis(time, shape=pars[1], scale=pars[2])
  } else if (dist=="lnorm") {
    stats::dlnorm(time, meanlog=pars[1], sdlog=pars[2])
  } else if (dist=="gamma") {
    stats::dgamma(time, shape=pars[1], rate=pars[2])
  } else if (dist=="gompertz") {
    flexsurv::dgompertz(time, shape=pars[1], rate=pars[2])
  } else if (dist=="gengamma") {
    flexsurv::dgengamma(time, mu=pars[1], sigma=pars[2], Q=pars[3])
  } else if (dist=="gengamma.orig") {
    flexsurv::dgengamma.orig(time, shape=pars[1], scale=pars[2], k=pars[3])
  } else {
    NA
  }
}

#' Calculate the value of the density function
#' @description Calculate the value of the density function, given specification as either parametric or Royston-Parmar splines model
#' @inheritParams calc_surv
#' @inherit calc_haz_par return
#' @export
#' @examples
#' calc_dens(time=1:5,
#'     type="spl",
#'     spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
#'     )
#' calc_dens(time=1:5,
#'     type="par",
#'     spec=list(dist="lnorm", pars=c(3,1))
#'     )
calc_dens <- function(time, type, spec){
  type <- tolower(substr(type, 1, 3))
  if (type=="spl") {
    flexsurv::dsurvspline(x = time,
                          gamma = spec$gamma,
                          beta = 0,
                          X = 0,
                          knots = spec$knots,
                          scale = spec$scale,
                          timescale = "log",
                          offset = 0)
  } else if (type=="par") {
    calc_dens_par(time, dist=spec$dist, pars=spec$pars)
  } else {
    NA
  }
}

#' Calculate restricted mean durations (parametric form)
#' @description Calculates the restricted mean duration, given a parametric statistical distribution.
#' @param dist is the statistical distribution (named per [flexsurv::flexsurvreg]).
#' @param pars is a vector of the parameters for that distribution.
#' - Exponential distribution (`exp`) requires the rate parameter.
#' - Weibull distribution (both `weibullPH` and `weibull` formulations) requires the shape and scale parameters.
#' - Log-logistic distribution (`llogis`) requires the shape and scale parameters.
#' - Log-normal distribution (`lnorm`) requires the meanlog and sdlog parameters.
#' - Gamma and Gompertz distributions (`gamma` and `gompertz`) require the shape and rate parameters.
#' - Generalized Gamma requires the mu, sigma and Q parameters if using the standard parameterization (`gengamma`) or shape, scale and k parameters if using the original parameterization (`gengamma.orig`).
#' @param Tw is the time period over which the restricted mean is calculated
#' @inherit calc_pdist seealso
#' @return the restricted mean duration, a numeric value.
#' @export
#' @examples
#' calc_rmd_par(20, "exp", 0.2)
#' calc_rmd_par(10, "lnorm", c(3, 1))
calc_rmd_par <- function(Tw, dist, pars) {
  if(dist=="exp") {
    flexsurv::rmst_exp(Tw, rate=pars[1], start=0)
  } else if (dist=="weibullPH") {
    flexsurv::rmst_weibullPH(Tw, shape=pars[1], scale=pars[2], start=0)
  } else if (dist=="weibull" | dist=="weibull.quiet") {
    flexsurv::rmst_weibull(Tw, shape=pars[1], scale=pars[2], start=0)
  } else if (dist=="llogis") {
    flexsurv::rmst_llogis(Tw, shape=pars[1], scale=pars[2], start=0)
  } else if (dist=="lnorm") {
    flexsurv::rmst_lnorm(Tw, meanlog=pars[1], sdlog=pars[2], start=0)
  } else if (dist=="gamma") {
    flexsurv::rmst_gamma(Tw, shape=pars[1], rate=pars[2], start=0)
  } else if (dist=="gompertz") {
    flexsurv::rmst_gompertz(Tw, shape=pars[1], rate=pars[2], start=0)
  } else if (dist=="gengamma") {
    flexsurv::rmst_gengamma(Tw, mu=pars[1], sigma=pars[2], Q=pars[3], start=0)
  } else if (dist=="gengamma.orig") {
    flexsurv::rmst_gengamma.orig(Tw, shape=pars[1], scale=pars[2], k=pars[3], start=0)
  } else {
    NA
  }
}

#' Calculate restricted mean durations
#' @description Calculates the restricted mean duration, given the form of a parametric distribution of Royston-Parmar splines
#' @param Tw is the time horizon (weeks) over which the mean should be calculated.
#' @param type is either "par" for regular parametric form (exponential, weibull etc) or "spl" for Royston-Parmar splines.
#' @param spec is a list comprising:
#' If type=="par": `dist` is the statistical distribution (named per [flexsurv::flexsurvreg]) and `pars` is a vector of the parameters for that distribution.
#' - Exponential distribution (`exp`) requires the rate parameter.
#' - Weibull distribution (both `weibullPH` and `weibull` formulations) requires the shape and scale parameters.
#' - Log-logistic distribution (`llogis`) requires the shape and scale parameters.
#' - Log-normal distribution (`lnorm`) requires the meanlog and sdlog parameters.
#' - Gamma and Gompertz distributions (`gamma` and `gompertz`) require the shape and rate parameters.
#' - Generalized Gamma requires the mu, sigma and Q parameters if using the standard parameterization (`gengamma`) or shape, scale and k parameters if using the original parameterization (`gengamma.orig`).
#' If type=="spl":
#' - `gamma` - Vector of parameters describing the baseline spline function, as described in [flexsurv::flexsurvspline]. This may be supplied as a vector with number of elements equal to the length of knots, in which case the parameters are common to all times. Alternatively a matrix may be supplied, with rows corresponding to different times, and columns corresponding to knots.
#' - `knots` - Vector of locations of knots on the axis of log time, supplied in increasing order. Unlike in [flexsurv::flexsurvspline], these include the two boundary knots.
#' - `scale` - Either "hazard", "odds", or "normal", as described in [flexsurv::flexsurvspline]. With the default of no knots in addition to the boundaries, this model reduces to the Weibull, log-logistic and log-normal respectively. The scale must be common to all times.
#' @param discrate Discounting rate (%) per year
#' @inherit calc_haz_par seealso
#' @inherit calc_haz_par return
#' @export
#' @examples
#' calc_rmd(Tw=200,
#'     type="spl",
#'     spec=list(gamma=c(0.1,0.2,0.1), knots=c(-5,2,4), scale="normal")
#'     )
#' calc_rmd(Tw=250,
#'     type="par",
#'     spec=list(dist="lnorm", pars=c(3,1))
#'     )
calc_rmd <- function(Tw, type, spec, discrate=0){
  vn <- (1+discrate)^(-convert_wks2yrs(Tw))
  type <- tolower(substr(type, 1, 3))
  if (type=="spl") {
    vn * flexsurv::rmst_survspline(t = Tw,
                                 gamma = spec$gamma,
                                 beta = 0,
                                 X = 0,
                                 knots = spec$knots,
                                 scale = spec$scale,
                                 timescale = "log",
                                 offset = 0
                                 )
  } else if (type=="par") {
    vn * calc_rmd_par(Tw=Tw, dist=spec$dist, pars=spec$pars)
  } else {
    NA
  }
}

#' Number of parameters used by parametric statistical distributions
#' @description Returns the number of parameters used by one or many statistical distributions, named as per flexsurv.
#' @param dist is the statistical distribution (named per [flexsurv::flexsurvreg])
#' @return a numeric value
#' @inherit calc_pdist seealso
#' @export
#' @examples
#' give_noparams_par("llogis")
give_noparams_par <- function(dist) {
  if (length(dist)!=1) stop("Multiple distributions entered, may call only one.")
  dplyr::case_when(
    dist=="exp" ~ 1,
    dist=="weibull" | dist=="weibullPH" | dist=="llogis" | dist=="lnorm" ~ 2,
    dist=="gamma" | dist=="gompertz" ~ 2,
    dist=="gengamma" | dist=="gengamma.orig" ~ 3,
    .default = NA
  )
}
  
#' Number of parameters used by parametric statistical distributions
#' @description Returns the number of parameters used by the specified statistical model
#' @inheritParams calc_surv
#' @return a numeric value
#' @inherit calc_pdist seealso
#' @export
#' @examples
#' give_noparams(type="par", spec=list(dist="weibullPH"))
#' give_noparams(type="spl", spec=list(gamma=c(1.1,2.1,3.1)))
give_noparams <- function(type, spec) {
  type <- tolower(substr(type, 1, 3))
  if (type=="spl") {
    length(spec$gamma)
  }
  else if (type=="par") {
    psm3mkv::give_noparams_par(dist=spec$dist)
  }
  else {
    NA
  }
}

#' Convert weeks to years
#' @description Convert weeks to years. Function allows hard-coding to be minimized elsewhere in the package.
#' @param weeks Number of weeks
#' @return Number of years = weeks x 7 / 365.25
#' @seealso [psm3mkv::convert_yrs2wks]
#' @export
convert_wks2yrs <- function(weeks) {
  weeks*7/365.25
}

#' Convert years to weeks
#' @description Convert weeks to years. Function allows hard-coding to be minimized elsewhere in the package.
#' @param years Number of years
#' @return Number of weeks = years / 7 * 365.25
#' @seealso [psm3mkv::convert_wks2yrs]
#' @export
convert_yrs2wks <- function(years) {
  years*365.25/7
}