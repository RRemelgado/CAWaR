#' @title meStack
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Stacking of raster layers with different extents
#' @param x A \emph{list} of \emph{RasterLayer} objects or a \emph{character} vector with the paths to \emph{raster} objects.
#' @param y A spatial object from which an extent can be derived.
#' @param z Object of class \emph{Date} with the acquisition date for each element in \emph{x}.
#' @param agg.fun Function used to aggregate images collected in the same date. Default is the mean.
#' @param derive.stats Logical argument. Default is FALSE.
#' @return A list containing a \emph{RasterStack} and related statistics.
#' @importFrom stats cor sd
#' @importFrom ggplot2 ggplot aes_string geom_ribbon geom_line
#' @importFrom lubridate is.Date
#' @importFrom raster raster extent extend res crs projectExtent projectRaster crop
#' @importFrom rgeos intersect
#' @details {The function stacks the raster objects specified in \emph{x}. For each element 
#' in \emph{x}, the function crops it by the extent of \emph{y} and, if their extents differ, 
#' fits the extent of \emph{x} to the one of \emph{y}. All new pixels are set to NA. If \emph{z} 
#' is provided, the function will then aggregate all bands acquired in the same date using the 
#' function provide with \emph{agg.fun}. If \emph{derive.stats} is set to TRUE, the function will 
#' return basic statistics for each band (i.e. min, max, mean and sd) together with a plot of the 
#' mean values. The final output of the function is a list containing:
#' \itemize{
#'  \item{\emph{stack} - \emph{RasterStack} object.}
#'  \item{\emph{dates} - Acquisition dates for each layer in \emph{stack}.}
#'  \item{\emph{image.stats} - Statistics for each band in the output \emph{RasterStack}.}
#'  \item{\emph{stats.plot} - Plot showing the mean, minimum and maximum values per band.}
#'  \item{\emph{control} - Logical vector showing which elements in \emph{x} where used to build the \emph{RasterStack}.}}}
#' @examples {
#' 
#' require(raster)
#' 
#' r1 <- raster(xmn=1, xmx=90, ymn=1, ymx=90, res=1, vals=1) # image 1
#' r2 <- raster(xmn=50, xmx=150, ymn=50, ymx=150, res=1, vals=1) # image 2
#' r0 <- raster(xmn=20, xmx=90, ymn=50, ymx=90, res=1, vals=1) # target extent
#' 
#' crs(r0) <- crs(r2) <- crs(r1)
#' 
#' mes <- meStack(list(r1, r2), r0)
#' plot(mes$stack)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

meStack <- function(x, y, z, agg.fun=mean, derive.stats=FALSE) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  # check if x is a valid object
  if (!class(x)[1] %in% c("character", "list")) {stop('"x" is not of a valid class')}
  if (length(x) == 1) {stop('"x" only has 1 element')}
  
  
  # read x as a list of rasters if it is a character vector with image paths
  if (is.character(x)) {
    x <- lapply(x, function(i) {
      r <- file.exists(i)
      if (r) {
        r <- try(raster(i), silent=TRUE)
        if (class(r)[1] != 'try-error') {return(r)} else {return(NULL)}
      } else {return(NULL)}})}
  
  # test if all rasters are usable
  if (sum(sapply(x, function(i) {is.null(i)})) > 0) {stop('one or more elements in "x" are not valid rasters')}
  
  # check if the pixel resolution of individual images differ
  pixel.res <- unique(sapply(x, function(i) {res(i)[1]}))
  if (length(pixel.res) > 1) {stop('the lements in "x" have more than 1 pixel resolution (same sensor?)')}
  
  # see if y is a valid spatial object
  e <- try(extent(y), silent=TRUE)
  if (class(e)[1] == 'try-error') {stop('"y" is not a valid spatial object')}
  y <- raster(extend(extent(y), pixel.res), res=pixel.res, crs=crs(y))

  # check if temporal information was provided
  if (!missing("z")) {
    if (!is.Date(z)) {stop('"z" is not a "Date" vector')}
    if (length(x) != length(z)) {stop('"x" and "z" have different lenghts')}
    tinfo <- TRUE
  } else {
    tinfo <- FALSE
  z <- vector("numeric", length(x))
  z[] <- NA}

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Build stack
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  tmp <- lapply(1:length(x), function(i) {
    
    ref <- projectExtent(y, crs(x[[i]])) # reference extent (used for reprojection)
    
    if (tryCatch(!is.null(intersect(x[[i]], ref)), error=function(e) return(FALSE), finally=TRUE)) {
        
      r <- projectRaster(crop(x[[i]], ref), y) # reprojects is necessary
      r <- extend(crop(r, y), y, value=NA) # crop/extend to extent
      ov <- list(image=r, date=z[i]) # return image/daye
      
    } else {ov <- NULL}
    
    return(ov)

  })

  c.var <- sapply(tmp, function(i) {return(!is.null(i))}) # control variable (which images where used?)
  o.stk <- stack(lapply(tmp[c.var], function(i) {i$image})) # output stack
  acqd <- do.call("c", lapply(tmp[c.var], function(i) {i$date})) # acquistion dates

  rm(tmp)
  
  # sort stack by time (if "z" is provided)
  if (tinfo) {

    ud <- unique(acqd) # unique dates

    # combine images with same dates
    o.stk <- stack(lapply(ud, function(d) {
      i <- which(acqd == d)
      if (length(i) > 1) {r <- calc(o.stk[[i]], agg.fun, na.rm=TRUE)} else {r <- o.stk[[i]]}
      return(r)}))

    si <- order(ud) # date sorting indices
    o.stk <- o.stk[[si]] # sort stack
    id <- ud[si] # sort dates

  } else {

    id <- 1:nlayers(o.stk) # image ID

  }

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Estimate statistics
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (derive.stats) {
    
    odf <- do.call("rbind", lapply(1:nlayers(o.stk), function(i) {
      
      s1 <- cellStats(o.stk[[i]], min, na.rm=TRUE)
      s2 <- cellStats(o.stk[[i]], max, na.rm=TRUE)
      s3 <- cellStats(o.stk[[i]], mean, na.rm=TRUE)
      s4 <- cellStats(o.stk[[i]], sd, na.rm=TRUE)
      
      return(data.frame(min=s1, max=s2, mean=s3, sd=s4))
      
    }))

    p <- ggplot(odf, aes_string(x="id")) + theme_bw() + geom_ribbon(aes_string(x='id', ymin='min', ymax='max'), fill="grey70") +
      geom_line(aes_string(y='mean')) + theme_bw() + xlab("\nBand ID") + ylab("Value\n")

  } else {

    odf <- NULL
    p <- NULL

  }

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. Derive output
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  return(list(stack=o.stk, dates=acqd, image.stats=odf, stats.plot=p, control=c.var))

}
