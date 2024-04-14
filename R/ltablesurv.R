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
# Functions relating to life table survival
# ltablesurv.R
# ==================================================================

#' VoneLOOKUP function
#' @description Function to lookup a single (one) value according to an index. Aims to behave similarly to VLOOKUP in Microsoft Excel.
#' @param oneindexval The single index value to be looked-up
#' @param indexvec The vector of indices to look-up in
#' @param valvec The vector of values corresponding to the vector of indices
#' @param method Method may be `floor`, `ceiling`, `arith` or `geom` (default).
#' @return Numeric value or vector, depending on the lookup/interpolation method chosen:
#' - `floor`: Floor value, where interpolation is required between measured values
#' - `ceiling`: Ceiling value, where interpolation is required between measured values
#' - `arith`: Arithmetic mean, where interpolation is required between measured values
#' - `geom`: Geometric mean, where interpolation is required between measured values
#' @seealso [psm3mkv::vlookup]
vonelookup <- function(oneindexval, indexvec, valvec, method="geom") {
  if (is.na(oneindexval)) return(NA) # Return NA rather than error if index value lookup is NA
  if (oneindexval<min(indexvec)) stop("Lookup value is below range of lookup table")
  if (oneindexval>max(indexvec)) stop("Lookup value is above range of lookup table")
  loc <- indexrange <- valrange <- NULL
  # Location of index values
  loc <- match(1, (oneindexval>=indexvec)*(oneindexval<=dplyr::lead(indexvec)))
  indexrange <- indexvec[loc:(loc+1)]
  valrange <- valvec[loc:(loc+1)]
  dplyr::case_when(
    method == "floor" ~ dplyr::if_else(oneindexval==indexrange[2], valrange[2], valrange[1]),
    method == "ceiling" ~ dplyr::if_else(oneindexval==indexrange[1], valrange[1], valrange[2]),
    method == "arith" ~ (valrange[1]*(indexrange[2]-oneindexval) + valrange[2]*(oneindexval-indexrange[1])) / (indexrange[2]-indexrange[1]),
    method == "geom" ~ ((valrange[1]^(indexrange[2]-oneindexval)) * (valrange[2]^(oneindexval-indexrange[1]))) ^ (1/(indexrange[2]-indexrange[1])),
    .default = ((valrange[1]^(indexrange[2]-oneindexval)) * (valrange[2]^(oneindexval-indexrange[1]))) ^ (1/(indexrange[2]-indexrange[1]))
  )
}

#' VLOOKUP function
#' @description Function to lookup values according to an index. Aims to behave similarly to VLOOKUP in Microsoft Excel.
#' @param indexval The index value to be looked-up (may be a vector of multiple values)
#' @param indexvec The vector of indices to look-up in
#' @param valvec The vector of values corresponding to the vector of indices
#' @param method Method may be `floor`, `ceiling`, `arith` or `geom` (default).
#' @return Numeric value or vector, depending on the lookup/interpolation method chosen:
#' - `floor`: Floor value, where interpolation is required between measured values
#' - `ceiling`: Ceiling value, where interpolation is required between measured values
#' - `arith`: Arithmetic mean, where interpolation is required between measured values
#' - `geom`: Geometric mean, where interpolation is required between measured values
#' @seealso [HMDHFDplus::readHMDweb] can be used to obtain lifetables from the Human Mortality Database
#' @export
#' @examples
#' # Suppose we have survival probabilities at times 0 to 20
#' times <- 0:20
#' survival <- 1-times*0.04
#' # We would like to look-up the survival probability at time 7
#' vlookup(7, times, survival)
#' # In this case, the floor, ceiling, arith and geom values are identical
#' # because survival time 7 is known, and no interpolation is necessary
#' vlookup(c(7, 7.5), times, survival)
#' # The second row of the returned tibble reveal different estimates of the survival at time 7.5.
#' # The values vary according to the interpolation method between
#' # observed survival values at times 7 and 8.
vlookup <- function(indexval, indexvec, valvec, method="geom") {
  # Checks that inputs are reasonable
  if (!is.numeric(indexvec)) stop("Index vector must be numeric")
  if (!is.numeric(valvec)) stop("Values vector must be numeric")
  if (length(indexvec)!=length(valvec)) stop("Index and values vectors must be same length")
  if (!all(indexvec==sort(indexvec))) stop("Index vector must sorted into ascending order")
  if (!all(indexvec==unique(indexvec))) stop("Index vector components must be unique")
  # Function for looking up value for one index value (oneindexval)
  sapply(indexval, vonelookup, indexvec=indexvec, valvec=valvec, method=method)
}

#' Calculate survival from a lifetable
#' @description Calculate survival from time zero to a given time, according to a provided lifetable
#' @param looktime The time(s) to which survival is to be estimated (from time zero).
#' @param lifetable The lifetable must be a dataframe with columns named `lttime` (years) and `lx`. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param method Method may be `floor`, `ceiling`, `arith` or `geom` (default).
#' @return Numeric survival probability
#' @export
#' @examples
#' ltable <- tibble::tibble(lttime=0:10, lx=10-(0:10))
#' calc_ltsurv(c(2, 2.5, 9.3), ltable)
calc_ltsurv <- function(looktime, lifetable=NA, method="geom"){
  if (!is.data.frame(lifetable)) {return(rep(1, length(looktime)))}
  if (lifetable$lttime[1]!=0) {stop("Lifetable must run from time zero")}
  vlookup(looktime, lifetable$lttime, lifetable$lx, method=method) / lifetable$lx[1]
}

#' Calculate mortality density from a lifetable
#' @description Calculate mortality density a given time, according to a provided lifetable
#' @param looktime The time(s) to which survival is to be estimated (from time zero).
#' @param lifetable The lifetable must be a dataframe with columns named `lttime` (years) and `lx`. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param method Method may be `floor`, `ceiling`, `arith` or `geom` (default). 
#' @return Numeric survival probability
#' @export
#' @examples
#' ltable <- tibble::tibble(lttime=0:10, lx=10-(0:10))
#' calc_ltdens(c(2, 2.5, 9.3), ltable)
calc_ltdens <- function(looktime, lifetable=NA, method="geom"){
  if (!is.data.frame(lifetable)) stop("Lifetable must be specified")
  if (lifetable$lttime[1]!=0) stop("Lifetable must run from time zero")
  # Floor time from lifetable
  tlo <- vlookup(looktime, lifetable$lttime, lifetable$lttime, method=method)
  pos <- match(tlo, lifetable$lttime)
  # Pick out useful lx values
  lx0 <- lifetable$lx[1]
  lxlo <- lifetable$lx[pos]
  lxhi <- lifetable$lx[pos+1]
  tdiff <- lifetable$lttime[pos+1]-lifetable$lttime[pos]
  # Average hazard
  hx <- -log(lxhi/lxlo)/tdiff
  # Density at looktime, f = h x S
  hx * lxlo / lx0
}

#' Calculate restricted life expectancy from a lifetable
#' @param Ty Time duration over which to calculate (default is 10 years). Assumes input is in years, and patient-level data is recorded in weeks.
#' @param lifetable The lifetable must be a dataframe with columns named `lttime` (years) and `lx`. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param discrate Discount rate (%) per year
#' @return List containing `ex_y` and `ex_w'`, the numeric (restricted) life expectancy in years and weeks respectively,
#' and `calcs`, a dataframe of the calculations.
#' @examples
#' # Create a lifetable. Must end with lx=0.
#' # ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
#' # calc_ex(lifetable=ltable, discrate=0.03)
#' # calc_ex(Ty=Inf, lifetable=ltable)
calc_ex <- function(Ty=10, lifetable, discrate=0) {
  lx <- lttime <- midtime <- vn <- beyond <- lxvn <- blxvn <- NULL
  # Lifetable must have minimum lx = 0
  if (min(lifetable$lx)!=0) stop("Lifetable must end with a zero lx")
  # Calculation
  res1 <- lifetable |>
    dplyr::mutate(
      midtime = dplyr::if_else(lx>0, (lttime + dplyr::lead(lttime))/2, 0),
      vn = (1+discrate)^(-midtime),
      lxvn = lx*vn,
      beyond = (lttime>Ty)*1,
      blxvn = beyond*lxvn
    )
  res2 <- res1 |>
    dplyr::summarize(
      lx0 = max(lxvn),
      Dx0 = sum(lxvn),
      Dxt = sum(blxvn)
    )
  ex <- (res2$Dx0-res2$Dxt)/res2$lx0
  list(
    ex_y = ex,
    ex_w = convert_yrs2wks(ex),
    calcs = res1)
}

#' Constrain survival probabilities according to hazards in a lifetable
#' Recalculated constrained survival probabilities (by week) as the lower of the original unadjusted survival probability and the survival implied by the given lifetable (assumed indexed as years).
#' @param survprob1 (Unconstrained) survival probability value or vector
#' @param survprob2 Optional survival probability value or vector to constrain on (default = NA)
#' @param lifetable Lifetable (default = NA)
#' @param timevec Vector of times corresponding with survival probabilities above
#' @return Vector of constrained survival probabilities
#' @export
#' @examples
#' ltable <- tibble::tibble(lttime=0:20, lx=c(1,0.08,0.05,0.03,0.01,rep(0,16)))
#' survprob <- c(1,0.5,0.4,0.2,0)
#' constrain_survprob(survprob, lifetable=ltable)
#' timevec <- 100*(0:4)
#' constrain_survprob(survprob, lifetable=ltable, timevec=timevec)
#' survprob2 <- c(1,0.45,0.35,0.15,0)
#' constrain_survprob(survprob, survprob2)
constrain_survprob <- function(survprob1, survprob2=NA, lifetable=NA, timevec=0:(length(survprob1)-1)) {
  # Check what exists
  s2exists <- !is.na(survprob2)[1]
  ltexists <- !is.na(lifetable)[1]
  # Survprob1 and survprob2 must have equal length (when survprob2 is defined)
  if(s2exists) {
    stopifnot("Survival probability vectors must have equal length" = length(survprob1)==length(survprob1))
    }
  # Vector of lifetables (or ones, if lifetable not specified)
  tN <- length(timevec)
  lxprob <- rep(1, tN)
  if(ltexists) {
    lxprob <- calc_ltsurv(convert_wks2yrs(timevec), lifetable)
  }
  # Cycle through each element
  adjsurv <- slx <- sprob <- rep(NA, tN)
  adjsurv[1] <- survprob1[1]
  for (i in 2:tN) {
    slx[i] <- ifelse(lxprob[i-1]==0, 1, lxprob[i]/lxprob[i-1])
    sprob[i] <- ifelse(survprob1[i-1]==0, 1, survprob1[i]/survprob1[i-1])
    sprob[i] <- ifelse(s2exists,
                       ifelse(survprob2[i-1]==0, sprob[i],
                                 pmin(sprob[i], survprob2[i]/survprob2[i-1])),
                        sprob[i])
    adjsurv[i] <- adjsurv[i-1] * pmin(slx[i], sprob[i])
  }
  return(adjsurv)
}