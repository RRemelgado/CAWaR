#' @title extractTS
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Extracts time series data from a \emph{RasterStack} for a \emph{SpatialPolygons} or a \emph{SpatialPolygonsDataFrame} object.
#' @param x Object of class \emph{SpatialPolygons}, \emph{SpatialPolygonsDataFrame}, \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param y A \emph{raster} object, a list of \emph{RasterLayer}'s or a numeric element.
#' @param z \emph{Numeric} vector with weights for each element in \emph{x} (when points).
#' @param id \emph{Numeric} vector with unique identifiers for \emph{x} (when points).
#' @return A \emph{list}.
#' @importFrom raster crs raster extent weighted.mean
#' @importFrom fieldRS poly2sample
#' @importFrom stats sd
#' @details {For each polygon in \emph{x} - if \emph{x} is a \emph{SpatialPolygons} and \emph{SpatialPolygonsDataFrame} object - the function 
#' identifies the overlapping pixels in \emph{y} and, for each pixel, estimates thepercentage area covered by the polygon. Using this data as 
#' weights, the function calculates the weighted mean for each band in \emph{y}. If \emph{y} is a numeric element, the function will build a 
#' raster with resolution equal to \emph{y} over which the pixel cover will be estimated. Moreover, if \emph{x} is a \emph{SpatialPoints} or 
#' a \emph{SpatialPointsDataFrame} object, the function will skip the pixel extraction step. In this case, the user may provide a vector with 
#' sample weights through \emph{z} and a vector of unique identifiers (reporting on e.g. the polygon membership) The function returns a list 
#' of three \emph{data.frame} objects where each row represents a different polygon in \emph{x}:
#' \itemize{
#' \item{\emph{pixel.info} - \emph{SpatialPointsDataFrame} with pixel-wise samples for each polygon (identified by the field \emph{id}).}
#' \item{\emph{polygon.info} - Mean, min, max and standard deviation of the pixel cover; centroid coordinates.}
#' \item{\emph{weighted.mean} - Weighted mean raster values (if \emph{y} is a raster object).}}}
#' @seealso \code{\link{analyzeTS}}
#' @examples {
#' 
#' require(raster)
#' require(fieldRS)
#' 
#' # read raster data
#' r <- brick(system.file("extdata", "ndvi.tif", package="fieldRS"))
#' 
#' # read field data
#' data(fieldData)
#' 
#' extractTS(fieldData[1:5,], r)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

extractTS <- function(x, y, z, id) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  # decides how to process x
  if(is.null(tryCatch(extent(x), error=function(e) return(NULL)))) {stop('"x" is not of a valid class')}
  if (class(x)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {is.pts <- FALSE}
  
  if (class(x)[1] %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
    if (missing(z)) {z <- rep(1, length(x))} else {
      if (length(z) != length(x)) {stop('"x" and "z" have different lengths')}
      if (!is.numeric(z)) {stop('"z" is not numeric')}}
    if (missing(id)) {id <- rep(1, length(x))} else {
      if (length(id) != length(x)) {stop('"x" and "id" have different lenghts')}
      if (!is.numeric(id)) {stop('"id" is not numeric')}}
    is.pts <- TRUE}
  
  # decides how to process y...
  if (!class(y) %in% c("numeric", "RasterStack", "RasterLayer", "RasterBrick", "list", "character")) {stop('"y" is not of a valid class')} else {
    
    # ... if it's numeric
    if (is.numeric(y)) {
      if (length(y) > 1) {stop('"y" has more than 1 element')}
      y <- raster(extent(x), res=y, crs=crs(x), vals=NA)}
    
    # ... if it's a raster
    if (class(y)[1] %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}
      nl <- nlayers(y)}
    
    # ... if it's a list
    if(is.list(y)) {
      if (sum(sapply(y, function(i) {is.raster(i)}))!=length(y)) {stop('one or more elements in "y" are not RasterLayers')}
      nl <- length(y)
    }
    
    # ... if it's a character vector
    if (is.character(y)) {
      y <- lapply(y, function(i) {tryCatch(raster(i), error=function(e) return(NULL))})
      if (sum(sapply(y, function(i) {is.null(i)})) > 0) {stop('one or more elements in "y" are not RasterLayers')}
    }
    
        
  }
    
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. derive samples
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (!is.pts) {
    
    # extract samples per polygon ID
    tmp <- lapply(1:length(x), function(i) {
      shp <- poly2sample(x[i,], y[[1]])
      return(list(shp=shp, id=replicate(length(shp), i)))})
    ind <- sapply(tmp, function(i) {return(!is.null(i$shp))})
    out.df <- do.call(rbind, lapply(tmp[ind], function(i) {i$shp@data}))
    out.df$id <- do.call("c", lapply(tmp[ind], function(i) {i$id}))
    out.df <- out.df[which(out.df$cover > 0),]
    
    x <- SpatialPointsDataFrame(out.df[,c("x","y")], out.df, proj4string=crs(x))
    
    rm(tmp, out.df)
    
  } else {
    
    x <- SpatialPointsDataFrame(x@coords, data.frame(x=x@coords[,1], y=x@coords[,2], cover=z, id=id), proj4string=crs(x))
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. function to estimate weighted mean (depends on data dimensions)
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # extract raster value
  if (!is.list(y)) {ev0 <- as.data.frame(extract(y, x))} else {ev0 <- as.data.frame(extract2(y, x))}
  
  # summary function
  summary.fun <- function(x, y) {apply(x, 2, function(j) {weighted.mean(j, y, na.rm=TRUE)})}

  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. estimate weighted mean
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # derive polygon statistics / temporal profiles
  out.val <- lapply(unique(x$id), function(i) {
    ind <- which(x$id == i)
    v <- summary.fun(ev0[ind,], x$cover[ind])
    odf <- data.frame(id=i, x=mean(x$x[ind]), y=mean(x$y[ind]), min.cover=min(x$cover[ind]), 
                      max.cover=max(x$cover[ind]), mean.cover=mean(x$cover[ind]), 
                      cover.sd=sd(x$cover[ind]), count=length(ind))
    return(list(val=v, info=odf))})
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 5. construct output
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # return list
  return(list(pixel.info=x, polygon.info=do.call(rbind, lapply(out.val, function(i) {i$info})), 
              weighted.mean=do.call(rbind, lapply(out.val, function(i) {i$val}))))

}
