#' @title matchIndices
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
#' x <- c(2200, 4500, 4600, 6400, 1600) # target
#' y <- c(1100, 1150, 1200, 6400, 1600) # reference
#' 
#' i <- matchIndices(x, y, 1) # find best match
#' 
#' # plot x (blue), and selected y (red)
#' plot(1:5, replicate(5,0), ylim=c(0,10000), type="l") 
#' lines(i$x, x[i$x], type="l", col="blue")
#' lines(i$y, y[i$y], type="l", col="red")
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

matchIndices <- function(x, y, z) {
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. Check variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!is.numeric(x)) {stop('"x" is not a numeric vector')}
  if (sum(!is.na(x)) == 0) {stop('no usable records in "x"')}
  if (!is.numeric(y)) {stop('"y" is not a numeric vector')}
  if (sum(!is.na(y)) == 0) {stop('no usable records in "y"')}
  if (!is.numeric(z)) {stop('"z" is not a numeric')}
  if (length(z) > 1) {stop('"z" has more than 1 element')}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Build distance matrix
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  d <- matrix(0, length(y), length(x))
  
  for (i in 1:length(y)) {
    
    for (j in 1:length(x)) {
      
      c <- c(d[(i-z),(j-z)], d[(i-z),j], d[i,(j-z)])
      if (sum(is.finite(c)) > 0) {d[i,j] <-  abs(y[i] - x[j]) + min(c[is.finite(c)])} else {d[i,j] <-  abs(y[i] - x[j])}
      
    }
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Find shortest path between x and y
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  cv <- vector("logical", length(x))
  oi <- vector("numeric", length(y))
  
  for (j in 1:length(y)) {
    
    i <- which(d[j,] == min(d[j,!cv], na.rm=TRUE) & !cv)[1]
    
    if (length(i) > 0) {
      
      oi[j] <- i
      cv[oi[j]] <- TRUE
      
    } else {
      
      oi[j] <- NA
      
    }
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. Return indices for non-NA pixels in x and y
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  i <- which(!is.na(oi))
  return(list(x=i, y=oi[i]))
  
}
