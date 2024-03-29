% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhoods.R
\name{convert_fit2spec}
\alias{convert_fit2spec}
\title{Obtain the type and specification as required in other package functions from a model fit}
\usage{
convert_fit2spec(fitsurv)
}
\arguments{
\item{fitsurv}{Fitted model (either parametric or Royston-Parmar splines model)}
}
\value{
List of type and specification
\itemize{
\item \code{type} is "spl" for splines model or "par" for parametric model
\item \code{spec} contains distribution (\code{dist}) and coefficients (\code{coefs}) if \code{type=="par"}
\item \code{spec} contains gamma values (\code{gamma}), knot locations (log scale, \code{knots}) and scale (\code{scale}) for Royston-Parmar splines model, if \code{type=="spl"}
}
}
\description{
Obtain the type and specification as required in other package functions from a model fit
}
\examples{
bosonc <- create_dummydata("flexbosms")
fits <- fit_ends_mods_spl(bosonc)
convert_fit2spec(fits$pfs[[3]]$result)
}
\seealso{
\code{\link[flexsurv:flexsurvspline]{flexsurv::flexsurvspline()}}
}
