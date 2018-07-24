#' @title phenoCropClass
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Spatialy explicit and phenology driven classification scheme for cropland mapping.
#' @param x A \emph{matrix} or \emph{data.frame}.
#' @param y A \emph{character} vector.
#' @param z A \emph{numeric} element. Default is 1.
#' @return A \emph{list} containing a set of reference profiles for each unique class in \emph{y}.
#' @importFrom stats cor
#' @details {Correlates \emph{x} with each row in \emph{y} using Dynamic Time Wraping (DTW) to match 
#' the time-series. \emph{z} sets the temporal buffer used to search to matching data points. The row 
#' in \emph{y} with the highest correlation is reported as the selected class. The final output is a 
#' \emph{data.frame} containing:
#' \itemize{
#'  \item{\emph{r2} - \eqn{R^{2}} between \emph{x} and each row \emph{y}.}
#'  \item{\emph{count} - Number of records used to estimate the \eqn{R^{2}}.}}}
#' @seealso \code{\link{analyzeTS}} \code{\link{phenoCropVal}}
#' @examples {
#' 
#' require(fieldRS)
#' 
#' # read reference profiles
#' data(referenceProfiles)
#' 
#' # target time series
#' x <- c(2200, 4500, 4600, 6400, 1600)
#' y <- referenceProfiles[,2:6]
#'
#' # Perform classification
#' c <-phenoCropClass(x, y)
#' head(c)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

phenoCropClass <- function(x, y, z) {
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. Check variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!is.numeric(x)) {stop('please provide "x" as a numeric vector')}
  if (!class(y)[1] %in% c("matrix", "data.frame")) {stop('"y" is not of a valid class')}
  if (length(x) != ncol(y)) {stop('"x" has a different lenght from the number of columns in "y"')}
  if (missing("z")) {z <- 1} else {if (!is.numeric(z)) {stop('"z" is not of a valid class')}}
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. Correlate time series
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # correlate x with reference profiles in y
  odf <- do.call(rbind, lapply(1:nrow(y), function(j) {
    
    i <- matchIndices(as.numeric(x), as.numeric(y[j,]), z)
    
    if (length(i) > 0) {
      
      r <- cor(as.numeric(x[i$x]), as.numeric(y[j,i$y]))^2
      i <- length(i$x)
      return(data.frame(r2=r, count=i))
      
    } else {return(data.frame(r2=NA, count=NA))}}))
  
  return(odf)
  
}
