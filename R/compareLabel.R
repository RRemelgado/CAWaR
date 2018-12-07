#' @title compareLabel
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Identifies samples that potentially require class label corrections.
#' @param x Object of class \emph{data.frame} with target profiles.
#' @param y Object of class \emph{data.frame} with reference profiles.
#' @param x.label \emph{Character vector} with classes of \emph{x}.
#' @param y.label \emph{Character vector} with classes of \emph{y}.
#' @param na.count Maximum number of NA values accepted.
#' @param max.length Maximum length of consecutive NA values accepted.
#' @return A \emph{list}.
#' @importFrom raster which.max
#' @importFrom stats ccf complete.cases
#' @importFrom stats na.omit
#' @details {The function cross-correlates \emph{x} and \emph{y}. Then,for each row, the function returns 
#' the element in \emph{y.label} with the highest correlation. The funal output of the function consists:
#' \itemize{
#'  \item{\emph{cross.cor} - Median, minimum and maximum values for each column in \emph{x} over each unique class in \emph{y}.}
#'  \item{\emph{label.compare} - \emph{data.frame} showing \emph{y.label} and the nbest match in \emph{y.label}.}
#'  \item{\emph{na.stats} - \emph{data.frame} showing the count and maximum number of consecutive NA values for each row in \emph{x}.}
#'  }
#'  Note that \emph{na.count} and \emph{max.length} determine which observations are judged. If These thresholds are exceeded, 
#'  the function will return NA.}
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
#' a.ts <- analyzeTS(as.data.frame(fieldDataTS$weighted.mean), fieldData$crop)
#' 
#' # extract reference profiles
#' rp <- as.data.frame(do.call(rbind, lapply(a.ts$y.statistics, function(i) {i$median})))
#' 
#' # compare labels
#' cl <- compareLabel(as.data.frame(fieldDataTS$weighted.mean), rp, fieldData$crop, a.ts$labels)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

compareLabel <- function(x, y, x.label, y.label, na.count=0, max.length=0) {
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (!class(x)[1] %in% c("matrix", "data.frame")) {stop('"x" is not of a valid class')}
  if (!class(y)[1] %in% c("matrix", "data.frame")) {stop('"y" is not of a valid class')}
  if (!is.character(x.label)) {stop('"x.label" is not a character vector')}
  if (!is.character(y.label)) {stop('"y.label" is not a character vector')}
  if (nrow(x) != length(x.label)) {stop('the number of rows in "x" is different from the length of "x.label"')}
  if (nrow(y) != length(y.label)) {stop('the number of rows in "y" is different from the length of "y.label"')}
  if (ncol(x) != ncol(y)) {stop('"x" and "y" hve a different number of rows')}
  if (!is.numeric(na.count)) {'"na.count" is not a numeric element'}
  if (!is.numeric(max.length)) {'"max.length" is not a numeric element'}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. evaluate distribution of NA values
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  na.stats <- do.call(rbind, lapply(1:nrow(x), function(i) {
    
    # ifnd segments of NA's
    seg <- rle(as.numeric(is.na(x[i,])))
    seg.id <- vector('numeric', length(seg$lengths))
    for (p in 1:length(seg$lengths)) {seg.id[(sum(seg$lengths[0:(p-1)])+1):sum(seg$lengths[1:p])] <- p}
    
    # evaluate segments
    uv <- unique(seg.id) # unique segments
    na.freq <- c(0, sapply(uv, function(u) {sum(seg.id==u)})) # check for the length of each segment
    
    # return the total number of NA's and the maximum NA segment size
    return(data.frame(na.count=sum(na.freq[seg$values==1]), max.length=max(na.freq[c(1, seg$values)==1])))
    
  }))
    
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. correlate target and reference time series
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  cc <- do.call(cbind, lapply(1:nrow(y), function(r) {
    
    return(apply(x, 1, function(i) {
      ind <- which(!is.na(i) & !is.na(y[r,]))
      if (length(ind) > 1) {
        return(cor(as.numeric(i[ind]),as.numeric(y[r, ind])))
      } else {return(NA)}}))
    
  }))
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. class matching
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  ind <- which(na.stats$na.count <= na.count & na.stats$max.length <= max.length)
  
  if (length(ind) == 0) {
    
    warning('no usable samples (too many NA values)')
    odf <- NULL
  
  } else {
    
    # find index for highest correlation
    tmp <- lapply(ind, function(i) {
      ii <- which.max(cc[i,])
      return(list(cl=y.label[ii], pc=cc[i,ii]))})
    
    class.label <- sapply(tmp, function(i) {i$cl}) # probable class
    pearson.cor <- sapply(tmp, function(i) {i$pc}) # probable class correlation 
    
    rm(tmp)  
    
    class.compare <- class.label == x.label # is the class correct?
    
    # build report data.frame
    odf <- data.frame(x.label=x.label, y.label=class.label, pearson=NA, compare=NA, stringsAsFactors=FALSE)
    odf$pearson[ind]=pearson.cor
    odf$compare[ind]=class.compare
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. report
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(cross.cor=cc, label.compare=odf, na.stats=na.stats))
  
}
