#' @title segmentRaster
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Segmentation of a raster object using k-means clustering and connected component labelling.
#' @param x Object of class \emph{RasterLayer}.
#' @param change.threshold Numeric element (0 - 100).
#' @param n.size Numeric element.
#' @return A \emph{list}.
#' @importFrom raster crs extent raster cellStats aggregate cellFromXY xyFromCell rowFromCell colFromCell getValues setValues rasterToPolygons
#' @importFrom stats kmeans
#' @details {The function clusters a raster object with k-means through \link[RStoolbox]{unsuperClass} and segments the output using
#' \code{\link{}ccLabel}. The number of samples used to determine the number of clusters can be defined through \emph{n.size} to reduce
#' the required computational time. The optimal number of clusters is determined with the elbow method. he elbow method looks at the
#' percentage of variance explained by the number of clusters. If adding a new cluster, refered here as k, leads to no improvement,
#' the function will use k-1 to derive. This breakpoint is determined by \emph{change.threshold} which controls the percent change
#' between the amount of variance explained by two consequent k values. the output is a list consisting of:
#'  \itemize{
#'  \item{\emph{class} - Cluster image.}
#'  \item{\emph{regions} - Segmented region image.}}}
#' @seealso \code{\link{}}
#' @examples {}
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

segmentRaster <- function(x, change.threshold=0, n.size=NULL) {

  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. Check input variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (!class(x)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {stop('"x" is not of a valid class')}
  if (!is.numeric(change.threshold)) {stop('"change.threshold" is not numeric')}
  if (length(change.threshold) > 1) {stop('"change.threshold" has more than 1 element')}
  if (change.threshold < 0 || change.threshold > 100) {stop('"change.threshold" is not within the valid range')}
  if (is.null(n.size)) {e <- getValues(x)} else {
    if (!is.numeric(n.size)) {stop('"n.size" should be numeric')}
    if (length(n.size) > 1) {stop('"n.size" has too many elements')}
    c <- calc(stack(lapply(1:nlayers(x), function(i) {return(!is.na(x[[i]]))})), sum)
    i <- which.max(c == nlayers(x))
    i <- sample(i, n.size)
    e <- extract(r.stk, i)
    rm(c, i)}

  e <- as.data.frame(e)
  e <- e[complete.cases(e),]

  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. cluster image and identify regions
  #-----------------------------------------------------------------------------------------------------------------------------------------------#

  # normalize image
  nb <- nlayers(x)
  for (b in 1:nb) {
    min.v <-cellStats(x[[b]], min, na.rm=TRUE) # minimum value
    max.v <- cellStats(x[[b]], max, na.rm=TRUE) # maximum value
    if (nb == 1) {x <- 2 * ((x[[b]] - min.v) / (max.v - min.v)) - 1} # normalize
    if (nb > 1) {x[[b]] <- 2 * ((x[[b]] - min.v) / (max.v - min.v)) - 1}
  }

  c <- 2 # number of clusters to split x into
  r <- 100 # control variable 1 (used in elbow method)
  d <- 100 # control variable 2 (used in elbow method)

  # find optimal number of clusters
  while (d > 0) {
    k = kmeans(e, c)
    k <- round((1-(k$betweenss / k$totss))*100)
    d <- r - k
    r <- k
    c <- c + 1}

  k1 <- unsuperClass(x, nSamples=n.size, nClasses=(c-1))$map # final cluster image
  k2 <- ccLabel(k1) # segmented image

  return(list(class=k1, regions=k2))

}
