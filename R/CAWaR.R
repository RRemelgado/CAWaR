#' CAWaR.
#'
#' @name CAWaR
#' @docType package
#' @import raster sp
NULL

#' Samples for fieldData
#'
#' Output of extractTS for fieldData.
#'
#' \itemize{
#'   \item{pixel.info}{Unique pixels covered by fieldData with information of x and y coordinates and cover percent.}
#'   \item{polygon.info}{Percent cover statistics and pixel count per polygon.}
#'   \item{weighted.mean}{Weighted-mean time-series for each polygon.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fieldDataTS
#' @usage data(fieldDataTS)
#' @format A data.frame
NULL

#' Region labels for fieldData
#'
#' Output of splitSamples for fieldData.
#'
#' \itemize{
#'   \item{region.id}{Region identifier showing which samples are grouped.}
#'   \item{region.frequency}{Frequency of samples per region.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fieldDataCluster
#' @usage data(fieldDataCluster)
#' @format A data.frame
NULL
