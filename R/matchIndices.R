#' @title matchIndices
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Matches two vectors with different lenghts based on their maximum value.
#' @param x Target numeric \emph{vector}.
#' @param y Reference numeric \emph{vector}.
#' @return A \emph{vector} with indices for \emph{x} of the same length as \emph{y}.
#' @details {The function starts by identifying the maximum values in \emph{x} and \emph{y}. Then, it counts the number 
#' of values before and after the maximum in \emph{y} and returns the indices for the same amount of values in \emph{y}.}
#' @examples {}
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

matchIndices <- function(x, y) {
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!is.numeric(x)) {stop('"x" is not a numeric vector')}
  if (!is.numeric(y)) {stop('"y" is not a numeric vector')}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Find indices
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # peak indices
  i0 <- which.max(x)
  i1 <- which.max(y)
  
  # ideal start and end indices for y
  i1s <- abs(i1-i0)+1
  i1e <- i1 + (length(x)-i0)
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 3. correct indices for x and y
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # start
  if (i1s < 1) {
    i1s <- 1
    i0s <- abs(i0-i1)+1
  } else {i0s <- 1}
  
  # end
  if (i1e > length(x)) {
    i1e <- length(x)
    i0e <- i0+(i1e-i1)
  } else {i0e <- length(y)}
  
  return(list(x=i0s:i0e, y=i0s:i0e))
  
}
