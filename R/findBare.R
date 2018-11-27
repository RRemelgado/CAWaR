#' @title findBare
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Determines if a vector is likely related to a 
#' @param x Numeric \emph{vector}.
#' @param y Change threshold (in %).
#' @importFrom raster which.max
#' @return A \emph{logical} value (TRUE if "bare" and FALSE if not).
#' @details {Checks if \emph{x} is a likely "bare" time-series. In other 
#' words, it evaluates if the time series if from a pixel that was likely 
#' not cultivated assuming that such pixels should have small seasonal 
#' increase in the NDVI relatively to cultivated pixels. The change threshold 
#' is defined by \emph{y}. The function starts by normalizing \emph{x} by 
#' itself avoiding case-specific thresholding. As a consequence, \emph{y} 
#' is specified as a percent value of change.}
#' @examples {
#' 
#' x <- c(293, 770, 1166, 1166, 1562, 2357, 3234, 
#' 5806, 5806, 5678, 5678, 5546, 5536, 5536, 5536, 
#' 5325, 5200, 4726, 3550, 2868, 2365, 2365, 2365)
#' 
#' n <- findBare(x, 100)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#


findBare <- function(x, y) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if (!is.numeric(x)) {stop('"x" is not a numeric vector')}
  x <- x[!is.na(x)] # remove missing values
  if (length(x) < 2) {stop('"x" has less than 2 missing values')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. Evaluate x
#---------------------------------------------------------------------------------------------------------------------#
  
  x <- x / max(x) # normalize x bz its maximum
  a <- max(x) - min(x) # amplitude of x
  return(a < y)
  
}
