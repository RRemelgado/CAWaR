#' @title extractTS
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Extracts time series data from a \emph{RasterStack} for a \emph{SpatialPoints} object.
#' @param x Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonsDataFrame}.
#' @param y A \emph{raster} object or a numeric element.
#' @return A \emph{list}.
#' @importFrom raster crs raster extent weighted.mean
#' @importFrom fieldRS poly2sample
#' @importFrom stats sd
#' @details {For each polygon in \emph{x}, the function identifies the overlapping pixels in \emph{y} and, for each pixel, estimates the
#' percentage area covered by the polygon. Using this data as weights, the function calculates the weighted mean for each band in \emph{y}. 
#' If \emph{y} is a numeric element, the function will build a raster with resolution equal to \emph{y} over which the pixel cover will be 
#' estimated. The function returns a list of three \emph{data.frame} objects where each row represents a different polygon in \emph{x}:
#' \itemize{
#' \item{\emph{pixel.info} - Extracted weighted-mean time-series.}
#' \item{\emph{polygon.info} - Mean, min, max and standard deviation of the pixel cover; centroid coordinates.}
#' \item{\emph{weighted.mean} - Weighted mean raster values (if \emph{y} is a raster object).}}}
#' @seealso \code{\link{analyzeTS}}
#' @examples {
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

extractTS <- function(x, y) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (!class(x)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {stop('"x" is not of a valid class')}
  
  if (!class(y) %in% c("numeric", "RasterStack", "RasterLayer", "RasterBrick")) {stop('"y" is not of a valid class')} else {
    
    if (class(y)[1] %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}
      nl <- nlayers(y)
      ev <- TRUE}
    
    if (is.numeric(y)) {
      if (length(y) > 1) {stop('"y" has more than 1 element')}
      y <- raster(extent(x), res=y, crs=crs(x))
      ev <- FALSE}
    
  }
    
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. derive samples
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  # extract samples per polygon ID
  tmp <- lapply(1:length(x), function(i) {
    shp <- poly2sample(x[i,], y[[1]])
    return(list(shp=shp, id=replicate(length(shp), i)))})
  ind <- sapply(tmp, function(i) {return(!is.null(i$shp))})
  out.df <- do.call(rbind, lapply(tmp[ind], function(i) {i$shp@data}))
  out.df$id <- do.call("c", lapply(tmp[ind], function(i) {i$id}))
  out.df <- out.df[which(out.df$cover > 0),]
  
  rm(tmp)
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. estimate weighted mean
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # if y is a raster object...
  if (ev) {
    
    # extract raster values
    ev0 <- as.data.frame(extract(y, out.df[,c("x", "y")]))
    
    # determine how values should be handled (different for single and multi-band rasters)
    if (nl >1) {
      summary.fun <- function(x,y) {apply(x, 2, function(j) {weighted.mean(j, y, na.rm=TRUE)})}
    } else {
      summary.fun <- function(x,y) {weighted.mean(x, y, na.rm=TRUE)}
    }
    
  }
  
  # derive polygon statistics / temporal profiles
  out.val <- lapply(unique(out.df$id), function(i) {
    ind <- which(out.df$id == i)
    if (ev) {v <- summary.fun(ev0[ind,], out.df$cover[ind])} else {v <- NULL}
    odf <- data.frame(id=i, x=mean(out.df$x[ind]), y=mean(out.df$y[ind]), min.cover=min(out.df$cover[ind]), 
                      max.cover=max(out.df$cover[ind]), mean.cover=mean(out.df$cover[ind]), 
                      cover.sd=sd(out.df$cover[ind]), count=length(ind))
    return(list(val=v, info=odf))})
  
  # return list
  return(list(pixel.info=out.df, polygon.info=do.call(rbind, lapply(out.val, function(i) {i$info})), 
              weighted.mean=do.call(rbind, lapply(out.val, function(i) {i$val}))))

}
