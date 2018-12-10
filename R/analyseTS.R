#' @title analyseTS
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Summarizes multi-band \emph{raster} data within each element of a \emph{SpatialPolygons} object.
#' @param x Object of class \emph{data.frame}.
#' @param y Vector of class \emph{character} or \emph{numeric} with a length equal to the number of rows in \emph{x}.
#' @param out.plot Specifies the data path where plots should be stored.
#' @return A \emph{SpatialPointsDataDrame} with the coordinate pairs for each of the sampeled pixels.
#' @importFrom ggplot2 ggplot aes_string theme_bw ggtitle geom_ribbon geom_line xlab ylab ylim ggsave
#' @importFrom stats median mad
#' @details {For each unique value in \emph{y}, the function will select the rows in \emph{x} that correspond to it and estimate the
#' median, Median Absolte Deviation (MAD), minimum, maximum, mean and standard deviation for each column. Then, the function will build
#' a plot showing the median and draw a buffer that expresses the minimum and maximum. The final output is a list consisting of:
#' \itemize{
#'  \item{\emph{y.statistics} - Median, minimum and maximum values for each column in \emph{x} over each unique class in \emph{y}.}
#'  \item{\emph{plots} - List of line plots for each unique element in \emph{y}.}}
#'  If \emph{out.plot} is set, the function will save each plot as 10x10 cm PNG files within the specified path.
#'  }
#' @seealso \code{\link{extractTS}} \code{\link{phenoCropVal}} \code{\link{phenoCropClass}}
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
#' data(fieldDataTS)
#' 
#' a.ts <- analyseTS(as.data.frame(fieldDataTS$weighted.mean), fieldData$crop)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

analyseTS <- function(x, y, out.plot=NULL) {
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!class(x)[1]%in%c('matrix', 'data.frame')) {stop('"x" is not of a valid class')}
  if (!class(y)[1]%in%c('numeric', 'character')) {stop('"y" is not of a valid class')}
  if (ncol(x) == 1) {warning('"x" has only 1 variable')}
  if (nrow(x) != length(y)) {stop('"x" and "y" have different lengths')}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. derive reference time series
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # base variables
  unique.y <- unique(y) # unique classes
  min.value <- min(apply(x, 1, min, na.rm=TRUE), na.rm=TRUE) # define data range (minimum)
  max.value <- max(apply(x, 1, max, na.rm=TRUE), na.rm=TRUE) # define data range (maximum)
  unique.id <- 1:ncol(x) # plot x value
  
  # derive statistics and plots
  tmp <- lapply(unique.y, function(u) {
    i <- which(y == u)
    v <- as.matrix(x[i,])
    d <- data.frame(v1=as.numeric(apply(x[i,], 2, median, na.rm=TRUE)), v2=as.numeric(apply(x[i,], 2, mad, na.rm=TRUE)),
                    v3=as.numeric(apply(x[i,], 2, min, na.rm=TRUE)), v4=as.numeric(apply(x[i,], 2, max, na.rm=TRUE)),
                    v5=as.numeric(apply(x[i,], 2, mean, na.rm=TRUE)), v6=as.numeric(apply(x[i,], 2, sd, na.rm=TRUE)), 
                    id=unique.id, class=u, stringsAsFactors=FALSE)
    colnames(d) <- c("median", "mad", "min", "max", "mean", "sd", "id", "class")
    d0 <- d[,c("median", "mad", "id")]
    d0$um <- d$median+d$mad
    d0$lm <- d$median-d$mad
    p <- ggplot(d0, aes_string(x="id")) + theme_bw() + ggtitle(u) + geom_ribbon(aes_string(x='id', ymin='lm', ymax='um'), fill="grey70") +
      geom_line(aes_string(y='median')) + theme_bw() + xlab("\nVariable ID") + ylab("Value\n") + ylim(min.value, max.value)
    return(list(stats=d, plot=p))})
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. return output as a list
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(labels=unique.y, y.statistics=lapply(tmp, function(i) {i$stats}), plots=lapply(tmp, function(i) {i$plot})))
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 5. save plots (if prompted)
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!is.null(out.plot)) {
    
    string.format <- paste0('%0', as.character(nchar(as.character(length(unique.y)))), 'd')
    
    for (p in 1:length(unique(unique.y))) {
      
      ggsave(file.path(out.plot, paste0(sprintf(string.format, p), '_', unique.y[p], '.png')), tmp[[p]]$plot, width=10, height=10, units="cm")
      
    }
    
  }
  
}
