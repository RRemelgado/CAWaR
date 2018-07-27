#' @title meStack
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Stacking of raster layers with different extents
#' @param x A \emph{list} of \emph{RasterLayers}.
#' @param y A spatial object from which an extent can be derived.
#' @param z Object of class \emph{Date} with the acquisition date for each element in \emph{x}.
#' @param mask Object of class \emph{RasterLayer} with the same extent as \emph{y}.
#' @param plot Logical argument.
#' @param fun A function. Default is mean.
#' @return A list containing a \emph{RasterStack} and related statistics.
#' @importFrom stats cor sd
#' @importFrom ggplot2 ggplot aes_string geom_ribbon geom_line
#' @details {The function stacks the raster objects specified in \emph{x}. For each element in \emph{x}, the function crops it by the
#' extent of \emph{y} and, if their extents differ, fits the extent of \emph{x} to the one of \emph{y}. All new pixels are set to NA. If
#' \emph{z} is provided, the function will then aggregate all bands acquired in the same date \emph{plot} is set to TRUE, the function will derive basic statistics for each band (i.e. min, max, mean, sd) as well as a plot showing
#' then mean, min and max for each band. If \emph{mask} is provided, the plot will be based on all non-NA pixels. The final output of the
#' function is a list containing:
#' \itemize{
#'  \item{\emph{stack} - \emph{RasterStack} object.}
#'  \item{\emph{statistics} - Statistics for each band in the output \emph{RasterStack}.}
#'  \item{\emph{plot} - Plot showing the mean, minimum and maximum values per band.}
#'  \item{\emph{control} - Logical vector showing which elements in \emph{x} where used to build the \emph{RasterStack}.}}}
#' @examples {
#' 
#' require(raster)
#' 
#' r1 <- raster(xmn=1, xmx=90, ymn=1, ymx=90, res=30)
#' r2 <- raster(xmn=50, xmx=150, ymn=50, ymx=150, res=30)
#' 
#' meStack(list(r1, r2), )
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

meStack <- function(x, y, z, mask=NULL, plot=FALSE, fun=mean) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (!is.character(x)) {stop('"x" is not a character vector')}

  e <- try(extent(y), silent=TRUE)
  if (class(e)[1] == 'try-error') {stop('"y" is not a valid spatial object')}

  if (exists("z")) {
    if (!is.Date(z)) {stop('"z" is not a "Date" vector')}
    if (length(x) != length(z)) {stop('"x" and "z" have different lenghts')}
    tinfo <- TRUE
  } else {
    tinfo <- FALSE
  z <- vector("numeric", length(x))
  z[] <- NA}

  if (!is.null(mask)) {if (extent(mask) != extent(y)) {stop('"y" and "mask" have different extents')}}

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Build stack
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  tmp <- lapply(1:length(x), function(i) {

    r <- raster(x[i]) # read data
    o <- checkOverlap(r, y) # check if usable
    if (o[1] > 0) {
      r <- extend(crop(r, y), y, value=NA) # crop/extend to extent
      c <- TRUE # control(used)
    } else {
      c <- FALSE
      image <- NULL
    } # control(not used)

    return(list(control=c, image=r, date=z[i]))

  })

  o.stk <- stack(lapply(tmp, function(i) {i$image})) # output stack
  c.var <- do.call("c", lapply(tmp, function(i) {i$control})) # control variable (which images where used?)
  acqd <- do.call("c", lapply(tmp, function(i) {i$date})) # acquistion dates

  # sort stack by time (if "z" is provided)
  if (tinfo) {

    ud <- unique(acqd) # unique dates

    # combine images with same dates
    o.stk <- stack(lapply(ud, function(d) {
      i <- which(acqd == d)
      if (length(i) > 1) {r <- calc(o.stk[[i]], fun, na.rm=TRUE)} else {r <- o.stk[[i]]}
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

  if (plot) {

    if (!is.null(mask)) {
      ind <- which.max(mask)
      sf <- function(x) {return(data.frame(mean=mean, min=min(x[ind], na.rm=TRUE),
                                           max=max(x[ind], na.rm=TRUE),
                                           sd=sd(x[ind], na.rm=TRUE)))}
    } else {
      sf <- function(x) {return(data.frame(mean=cellStats(x, mean, na.rm=TRUE),
                                           min=cellStats(x, min, na.rm=TRUE),
                                           max=cellStats(x, max, na.rm=TRUE),
                                           sd=cellStats(x, sd, na.rm=TRUE)))}
    }

    odf <- do.call("rbind", lapply(1:nlayers(o.stk), function(i) {sf(o.stk[[i]])}))
    odf$id <- id

    p <- ggplot(mv, aes_string(x="id")) + theme_bw() + geom_ribbon(aes_string(x='id', ymin='min', ymax='max'), fill="grey70") +
      geom_line(aes_string(y='mean')) + theme_bw() + xlab("\nBand ID") + ylab("Value\n")

  } else {

    odf <- NULL
    p <- NULL

  }

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. Derive output
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  return(list(stack=o.stk, statistics=odf, plot=p, id=id))

}
