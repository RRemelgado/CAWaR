#' @title compareLabel
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Identifies samples that potentially require class label corrections.
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
#'  \item{\emph{na.stats} - \emph{data.frame} showing the count and maximum segment size of NA values for each row in \emph{x}.}
#'  }}
#' @seealso \code{\link{extractTS}} \code{\link{phenoCropVal}} \code{\link{phenoCropClass}}
#' @examples {
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

compareLabel <- function(x, y, x.label, y.label, na.count=0, max.na.count=0) {
  
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
# 2. evaluate distribution of NA values
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  na.stas <- do.call(rbind, lapply(1:nrow(x), function(i) {
    
    # ifnd segments of NA's
    seg <- rle(is.na(x))
    seg.id <- vector('numeric', length(pd))
    for (p in 1:length(seg$lengths)) {seg.id[(sum(seg$lengths[0:(p-1)])+1):sum(seg$lengths[1:p])] <- p}
    
    # evaluate segments
    uv <- unique(seg.id) # unique segments
    na.freq <- sapply(uv, function(u) {sum(seg.id==u)}) # check for the length of each segment
    
    # return the total number of NA's and the maximum NA segment size
    return(data.frame(na.count=sum(count.na[seg$values==1]), na.max.count=max(count.na[seg$values==1])))
    
  }))
    
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. correlate target and reference time series
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  cc <- ccf(x, y, type="cross-correlation", na.action=na.omit)
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. class matching
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  
  ind <- which(na.statss$na.count <= na.count & na.stats$na.max.count <= na.max.count)
  
  if (length(ind) == 0) {
    
    warning('no usable samples (too many NA values)')
    odf <- NULL
  
  } else {
    
    # find index for highest correlation
    tmp <- lapply(ind, function(i) {
      ind <- which.max(cc[i,])
      return(cl=y.label[ind], pc=cc[i,ind])})
    
    class.label <- sapply(tmp, function(i) {i$cl}) # probable class
    pearson.cor <- sapply(tmp, function(i) {i$pc}) # probable class correlation 
    
    rm(tmp)  
    
    class.compare <- class.label == x.label # is the class correct?
    
    # build report data.frame
    odf <- data.frame(x=x.label, y=class.label, peason=NA, compare=NA, stringsAsFactors=FALSE)
    odf$pearson[ind]=pearson.cor
    odf$compare[ind]=class.compare
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. report
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(cross.cor=cc, label.compare=odf, na.stats=na.stats))
  
}