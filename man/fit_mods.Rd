% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_mods}
\alias{fit_mods}
\title{Fit survival regressions with multiple splines or parametric distributions}
\usage{
fit_mods(durn1, durn2 = NA, evflag, type, spec)
}
\arguments{
\item{durn1}{Start time}

\item{durn2}{End time}

\item{evflag}{Event flag}

\item{type}{Type of model ("spl" for spline or "par" for parametric)}

\item{spec}{Specification of model}
}
\value{
A list by distribution, each containing two components:
\itemize{
\item result: A list containing information about the fitted model.
\item error: Any error message returned on fitting the regression (NULL indicates no error).
}
}
\description{
Fit survival regressions with multiple splines or parametric distributions
}
