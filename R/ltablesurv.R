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

# Mortality data - England & Wales, Female, 2019
# table <- HMDHFDplus::readHMDweb(CNTRY="GBRTENW",
#                                 item="fltper_1x1",
#                                 username="dom.muston@gmail.com",
#                                 password="Z8W@DJyMqv"
#) |>
#  dplyr::filter(Year==2019) |>
#  dplyr::select(Age, lx) |>
#  dplyr::mutate(Timew = ceiling((Age-50)*365.25/7)) |>
#  dplyr::filter(Timew>=0)|>
#  dplyr::select(-Age)
# ltable <- tibble::tibble(Timew=mort$Timew, lx=mort$lx)
# calc_ltsurv(0:10, ltable)

#' VLOOKUP function
#' @description Function to lookup values according to an index. Aims to behave similarly to VLOOKUP in Microsoft Excel.
#' @param indexval The index value to be looked-up (may be a vector of multiple values)
#' @param indexvec The vector of indices to look-up in
#' @param valvec The vector of values corresponding to the vector of indices
#' @return A tibble where each row provides results for each index value. Columns are:
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
#' # The values vary according to the interpolation method between observed survival values at times 7 and 8.
vlookup <- function(indexval, indexvec, valvec) {
  # Checks that inputs are reasonable
  if (!is.numeric(indexvec)) stop("Index vector must be numeric")
  if (!is.numeric(valvec)) stop("Values vector must be numeric")
  if (length(indexvec)!=length(valvec)) stop("Index and values vectors must be same length")
  if (!all(indexvec==sort(indexvec))) stop("Index vector must sorted into ascending order")
  if (!all(indexvec==unique(indexvec))) stop("Index vector components must be unique")
  ret <- NULL
  # Function for looking up value for one index value (oneindexval)
  onelookup <- function(oneindexval) {
    stopifnot(oneindexval >= min(indexvec), oneindexval<=max(indexvec))
    loc <- indexrange <- valrange <- NULL
    # Location of index values
    loc <- match(1, (oneindexval>=indexvec)*(oneindexval<=dplyr::lead(indexvec)))
    indexrange <- indexvec[loc:(loc+1)]
    valrange <- valvec[loc:(loc+1)]
    list(
      floor = dplyr::if_else(oneindexval==indexrange[2], valrange[2], valrange[1]),
      ceiling = dplyr::if_else(oneindexval==indexrange[1], valrange[1], valrange[2]),
      arith = (valrange[1]*(indexrange[2]-oneindexval) + valrange[2]*(oneindexval-indexrange[1])) / (indexrange[2]-indexrange[1]),
      geom = ((valrange[1]^(indexrange[2]-oneindexval)) * (valrange[2]^(oneindexval-indexrange[1]))) ^ (1/(indexrange[2]-indexrange[1]))
    )
  }
  ret <- matrix(unlist(sapply(indexval, onelookup)), ncol=length(indexval))
  ret <- tibble::tibble(floor=ret[1,], ceiling=ret[2,], arith=ret[3,], geom=ret[4,])
  # Location
  return(ret)
}

#' Calculate survival from a lifetable
#' @description Calculate survival from time zero to a given time, according to a provided lifetable
#' @param time The time(s) to which survival is to be estimated (from time zero).
#' @param lifetable The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @return Numeric survival probability
#' @export
#' @examples
#' ltable <- tibble::tibble(time=0:10, lx=1-time*0.08)
#' calc_ltsurv(c(2, 2.5, 11), ltable)
calc_ltsurv <- function(time, lifetable){
  if (lifetable$time[1]!=0) stop("Lifetable must run from time zero")
  vlookup(time, lifetable$time, lifetable$lx)$geom / lifetable$lx[1]
}