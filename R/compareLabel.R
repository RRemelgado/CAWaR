#' @title compareLabel
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Identifies samples tha tpotentially need class label corrections.
#' @param x Object of class \emph{data.frame} with target profiles.
#' @param y Object of class \emph{data.frame} with reference profiles.
#' @param x.lab \emph{Character vector} with classes of \emph{x}.
#' @param y.lab \emph{Character vector} with classes of \emph{y}.
#' @return A \emph{list}.
#' @importFrom raster which.max
#' @importFrom stats ccf complete.cases
#' @details {The function cross-correlates \emph{x} and \emph{y}. Then,for each row, the function returns 
#' the element in \emph{y.label} with the highest correlation. The funal output of the function consists:
#' \itemize{
#'  \item{\emph{cross.cor} - Median, minimum and maximum values for each column in \emph{x} over each unique class in \emph{y}.}
#'  \item{\emph{label.comapre} - \emph{data.frame} showing \emph{y.label} and the nbest match in \emph{y.label}.}
#'  }
#' @seealso \code{\link{extractTS}} \code{\link{phenoCropVal}} \code{\link{phenoCropClass}}
#' @examples {
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

compareLabel <- profile(x,y, x.label, y.label) {
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!class(x)[1] %in% c("matrix", "data.frame")) {stop('"x" is not of a valid class')}
  if (!class(y)[1] %in% c("matrix", "data.frame")) {stop('"y" is not of a valid class')}
  if (!is.character(x.label)) {stop('"x.label" is not a character vector')}
  if (!is.character(y.label)) {stop('"y.label" is not a character vector')}
  if (norw(x) != length(x.label)) {stop('the number of rows in "x" is different from the length of "x.label"')}
  if (norw(y) != length(y.label)) {stop('the number of rows in "y" is different from the length of "y.label"')}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. correlate target and reference time series
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  cc <- ccf(x, y, type="cross-correlation")
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. class matching
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  # find index for highest correlation
  tmp <- lapply(1:nrow(x), function(i) {
    ind <- which.max(cc[i,])
    return(cl=y.label[ind], pc=cc[i,ind])})
  
  class.label <- sapply(tmp, function(i) {i$cl}) # probable class
  pearson.cor <- sapply(tmp, function(i) {i$pc}) # probable class correlation
  
  rm(tmp)  
    
  class.compare <- class.label == x.label # is the class correct?
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. report
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  odf <- data.frame(x=x.label, y=class.label, pearson=pearson.cor, compare=class.compare, stringsAsFactors=FALSE)
  
  return(list(cross.cor=cc, label.compare=odf))
  
}