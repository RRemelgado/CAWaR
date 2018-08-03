#' @title countCropCycles
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Matches two vectors with different lenghts based on their maximum value.
#' @param x Target numeric \emph{vector}.
#' @param y Reference numeric \emph{vector}.
#' @param z A \emph{numeric} element.
#' @importFrom raster which.max
#' @return A \emph{list} with selected indices for \emph{x} and \emph{y}. 
#' @details {Uses Dynamic Time Wrapping (DTW) to match \emph{x} and \emph{y}. \emph{z} determines 
#' the buffer size - expressed in number of data points - used to search for matching records.}
#' @examples {
#' 
#' x <- c(293, 770, 1166, 1166, 1562, 2357, 3234, 
#' 5806, 5806, 5678, 5678, 5546, 5536, 5536, 5536, 
#' 5325, 5200, 4726, 3550, 2868, 2365, 2365, 2365)
#' 
#' n <- ountCropCycles(x)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

countCropCycles <- function(x) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if(!is.numeric(x)) {stop('"x" is not numeric')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. classify crop cylces
#---------------------------------------------------------------------------------------------------------------------#
  
  x1 <- x-mean(x) # identify break-point
  x1[x1 > 0] <- 1 # points above break-point (crop maturity)
  x1[x1 < 0] <- 0 # points below break-point (recently cultivated)
  s = rle(as.numeric(x1)) # identify growth cycles
  n <- sum(s$values == 1) # count cycles
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. identify breakpoints
#---------------------------------------------------------------------------------------------------------------------#
  
  s.id <- vector("numeric", length(x))
  for (i in 1:length(s$lengths)) {s.id[(sum(s$lengths[0:(i-1)])+1):sum(s$length[1:i])] <- i}
  
  
  return(n)
  
}
