#' @title splitSamples
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Aggregates a spatial object into regions.
#' @param x A \emph{SpatialPoints} or a \emph{SpatialPolygons} object.
#' @param y A \emph{RasterLayer}.
#' @param z A vector.
#' @param agg.radius Numeric element.
#' @return A list.
#' @importFrom raster res crs rasterize which.max rowFromCell colFromCell extract cellStats
#' @importFrom rsMove checkOverlap
#' @importFrom fieldRS spCentroid ccLabel
#' @details {For each class in \emph{z}, the function converts the elements in \emph{x} into a raster layer using \emph{y} as a basis. Then, 
#' it aggregates all pixels that are within a given distance of each other - defined by \emph{agg.radius} using \code{\link{ccLabel}}. The 
#' output is a list consisting of:
#' \itemize{
#'  \item{\emph{region.id} - Class dependent region label for each element in \emph{x}.}
#'  \item{\emph{region.frequency} - Pixel count for each unique value in \emph{region.id}.}}}
#' @seealso \code{\link{phenoCropVal}} \code{\link{phenoCropClass}}
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
#' fieldData <- fieldData[3:4,]
#' 
#' # find polygon clusters
#' k <- splitSamples(fieldData, r, fieldData$crop, agg.radius=30)
#' fieldData$ID <- as.factor(k$region.id)
#' 
#' # plot regions with labels
#' spplot(fieldData["ID"])
#' 
#' # show pixel count per region
#' head(k$region.frequency)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

splitSamples <- function(x, y, z, agg.radius=agg.radius) {
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. Check input variables 
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!class(x)[1] %in% c('SpatialPoints', 'SpatialPointsDataFrame', 'SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('"x" is not a valid spatial object')}
  if (!class(y) %in% c('RasterLayer', 'RasterStack', 'RasterBrick')) {stop('"y" is not a valid raster object')}
  
  if (class(x)[1] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {centroids <- spCentroid(x)} else {centroids <- x@coords}
  
  if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}
  if (checkOverlap(x,y)[1] != 100) {stop('"x" is not contained by "y"')}
  rdims <- dim(y) # raster dimensions
  
  if (exists("z")) {if (!is.vector(z)) {stop('"z" should be a vector')}} else {z <- replicate(length(x), 1)}
  unique.z <- unique(z) # target classes
  
  if (!is.numeric(agg.radius)) {stop('"agg.radius" is not a numeric element')}
  if (length(agg.radius) > 1) {stop('"agg.radius" has more than 1 element')}
  agg.radius <- round(agg.radius/res(y)[1]) # number of pixels
  if (agg.radius == (round(33 / 2)*2)) {agg.radius <- agg.radius + 1}
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. derive region indices
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  region.id <- vector('character', length(z)) # region id for each sample
  region.freq <- vector('list', length(z)) # region pixel count
  
  for (c in 1:length(unique.z)) {
    
    ri <- which(z == unique.z[c]) # target samples for class
    regions <- extend(rasterize(x[ri,], crop(y[[1]],x[ri,]), field=1, background=NA), agg.radius) # sample mask
    regions <- focal(regions, matrix(0,agg.radius, agg.radius), function(j) {sum(!is.na(j))}) > 0 # dilate
    regions <- crop(regions, y)
    
    # label regions
    regions <- ccLabel(regions)$regions
    
    # update region id's
    region.id[ri] <- paste0(unique.z[c], "_", sprintf("%003d", extract(regions, centroids[ri,1:2])))
    
    # count pixels per region
    urv <- unique(regions)
    urv <- urv[urv > 0]
    region.freq[[c]] <- data.frame(id=urv, count=sapply(urv, function(r) {cellStats(regions==r, sum)}))
    region.freq[[c]]$id <- paste0(unique.z[c], "_", sprintf("%003d", region.freq[[c]]$id))
    
  }
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 3. return output
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(region.id=region.id, region.frequency=do.call(rbind, region.freq)))
  
}
