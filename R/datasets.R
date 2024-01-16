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
# ======================================================================
#
# Datasets.R
# These functions are used to create dummy datasets to illustrate package use
# create_dummydata
# create_dummydata_survcan
# create_dummydata_flexbosms
#
# ======================================

#' Create dummy dataset for illustration
#' @description Create dummy dataset to illustrate [psm3mkv]
#' @param dsname Dataset name, as follows:
#' * 'flexbosms' provides a dataset based on [flexsurv::bosms3()]. This contains all the fields necessary for [psm3mkv].
#' * 'survcan' provides a dataset based on [survival::cancer()]. This contains the necessary ID and overall survival fields. You will additionally need to supply PFS and TTP data (fields pfs.durn, pfs.flag, ttp.durn and ttp.flag) to use [psm3mkv].
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @export
#' @examples
#' create_dummydata("survcan") |> head()
#' create_dummydata("flexbosms") |> head()
create_dummydata <- function(dsname) {
  if (dsname=="survcan") {create_dummydata_survcan()}
  else if (dsname=="flexbosms") {create_dummydata_flexbosms()}
  else {stop("Incorrect dataset specified. Must be survcan or flexbosms.")}
}

#' Create survcan dummy dataset for illustration
#' @description Create 'survcan' dummy dataset to illustrate [psm3mkv]
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @seealso [create_dummydata()]
#' @importFrom rlang .data
create_dummydata_survcan <- function() {
  # Declare local variables
  ptid <- os.durn <- os.flag <- NULL
  # Create dataset
  survival::cancer |>
    dplyr::mutate(
      ptid = dplyr::row_number(),
      os.durn = .data$time/7,
      os.flag = .data$status-1
    ) |>
    dplyr::select(ptid, os.durn, os.flag)
}

#' Create flexbosms dataset for illustration
#' @description Create 'flexbosms' dummy dataset to illustrate [psm3mkv]
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @seealso [create_dummydata()]
#' @importFrom rlang .data
create_dummydata_flexbosms <- function() {
  # Declare local variables
  id <- pfs.durn <- pfs.flag <- NULL
  os.durn <- os.flag <- ttp.durn <- ttp.flag <- NULL
  # Create dataset
  flexsurv::bosms3 |>
        tidyr::pivot_wider(id_cols = "id",
                     names_from = "trans",
                     values_from = c("Tstart", "Tstop", "status")
        ) |>
        dplyr::mutate(
          conv = 365.25/(12*7), # Convert from months to weeks
          pfs.durn = .data$conv * .data$Tstop_2,
          pfs.flag = .data$status_1 + .data$status_2,
          os.durn = .data$conv * ifelse(is.na(.data$Tstop_3), .data$Tstop_2, .data$Tstop_3),
          os.flag = ifelse(is.na(.data$Tstop_3), .data$status_2, .data$status_3),
          ttp.durn = .data$pfs.durn,
          ttp.flag = .data$status_1
        ) |>
        dplyr::select(id, pfs.durn, pfs.flag,
                      os.durn, os.flag, ttp.durn, ttp.flag) |>
        dplyr::rename(ptid = "id")
}
