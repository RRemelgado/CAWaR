#' @title findBare
#---------------------------------------------------------------------------------------------------------------------#
#' @description Determines if a vector is likely related to unused non-croped areas. 
#' @param x Numeric \emph{vector}.
#' @param y Change threshold (in %)
#' @return A logical element.
#' @details {Judges if \emph{x} is a time series of non-agricultural land (poitnetial bare land) based on its amplitude.}
#' @export

#---------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------#

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
