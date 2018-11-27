#' @title extract2
#----------------------------------------------------------------------------------------------------------------------------------#
#' @description Extract of values from multi-extent raster objects with a spatial object. 
#' @param x A \emph{character} vector with the paths to \emph{RasterLayers} or a \emph{list} of \emph{RasterLayers}.
#' @param y An object of class \emph{SpatialPoints} or \emph{SpatialPolygons}.
#' @param x.date Object of class \emph{Date} with the acquisition dates of each element in \emph{x}.
#' @param out.date Object of class \emph{Date} with the desired output dates.
#' @param time.buffer Two-element, numeric vector. 
#' @return A \emph{list} object.
#' @importFrom raster raster
#' @importFrom lubridate is.Date
#' @importFrom rsMove intime
#' @details {Creates a rectangular fishnet in a \emph{SpatialPolygon} format based on the 
#' extent of \emph{x} and the value of \emph{y} which defines the spatial resolution.}
#' @export

#----------------------------------------------------------------------------------------------------------------------------------#

extract2 <- function(x, y, x.date, out.date, time.buffer=c(365,365)) {
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#----------------------------------------------------------------------------------------------------------------------------------#
  
  # is x composed of readble raster data?
  if (!is.character(x) | !is.list(x)) {stop('"x" is not of a valid class')}
  if (is.character(x)) {x <- lapply(x, function(i) {return(tryCatch(raster(i), error=function(e) return(NULL)))})}
  if (is.list(x)) {x <- lapply(x, function(i) {return(tryCatch(is.raster(i), error=function(e) return(NULL)))})}
  if (sum(sapply(x, function(i) {is.null(i)) > 0})) {stop('one or more elements in "x" are not valid (check entries)')}
  
  # is y a shapefile?
  if (!tryCatch(extent(y), error=function(e) return(FALSE))) {stop('"y" is not a valid spatial object')}
  
  # is time information provided? (required for interpolation)
  if (!missing(x.date)) {
    
    # date information
    if (!is.Date(x.date)) {stop('"x.date" provided but not a Date object'}
    if (length(x.date) != length(x)) {stop('"x.date" and "x" should have the same lenght')}
    if (!missing(out.date)) {
      if (!is.Date(out.date)) {stop('"out.date" provided but not a Date object')}
      c1 <- sum(x.date >= min(out.date) & x.date <= max(out.date))
      c2 <- sum(out.date >= min(x.date) & out.date <= max(x.date))
      if (c1 == 0 & c2 == 0) {stop('"out.date" is not contained by "x.date (cannot interpolate)')}
    } else {out.date <- x.date}
    int.timne <- TRUE
    
    # temporal buffer to search for dates to interpolate from
    if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}
    if (length(time.buffer)==1) {time.buffer <- c(time.buffer, time.buffer)}
    
  } else {int.time <- FALSE}
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 2. extract raster values  
#----------------------------------------------------------------------------------------------------------------------------------#
  
  out.df <- do.call(cbind, lapply(x, function(i) {return(extract(i, y))}))
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 3. interpolate missing values (when prompted)  
#----------------------------------------------------------------------------------------------------------------------------------#
  
  if (int.time) {out.df <- intime(as.matrix(out.df), as.numeric(x.dates), as.numeric(out.dates), time.buffer)}
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 4. derive output
#----------------------------------------------------------------------------------------------------------------------------------#

  return(list(values=out.df, dates=out.dates))

}