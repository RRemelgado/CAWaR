#' @title findBare
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Shows is a vector is related to bare/fallow land..
#' @param x A numeric \emph{vector}.
#' @param y A numeric element.
#' @return A logical element.
#' @details {The function computes the amplitude of \emph{x} and returns TRUE if this value is below the threshold given by \emph{y}.}
#' @examples {
#' 
#' # vectors to test
#' x1 <- c(0.1,0.2,0.1,0.2,0.4,0.6,0.8,0.4,0.1) # simulated, crop profile
#' x2 <- c(0.1,0.3,0.2,0.2,0.1,0.3,0.1,0.2,0.2) # simulated, non-crop profile
#' 
#' # compare profiles
#' plot(x1, type="l")
#' lines(x2, col="red")
#' 
#' findBare(x1, 100) # returns TRUE
#' findBare(x2, 100) # returns FALSE
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

findBare <- function(x, y) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input
#---------------------------------------------------------------------------------------------------------------------#
  
  if (!is.numeric(x)) {stop('"x" is not a numeric vector')}
  x <- x[!is.na(x)] # remove missing values
  if (length(x) < 2) {stop('"x" has less than 2 missing values')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. Evaluate x
#---------------------------------------------------------------------------------------------------------------------#
  
  return((max(x) - min(x)) < y)
  
}
  