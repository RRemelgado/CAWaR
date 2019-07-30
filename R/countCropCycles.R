#' @title countCropCycles
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Matches two vectors with different lengths based on their maximum value.
#' @param x Target numeric \emph{vector}.
#' @param min.length Two element numeric \emph{vector}.
#' @importFrom raster which.max
#' @return A \emph{numeric} element with the number of crop cycles in \emph{x}. 
#' @details {The function counts the number of value segments in \emph{x} that are above its mean 
#' effectively counting the number of crop cycles. Before reporting the final value, \emph{min.length} 
#' is used to filter outliers. The first element filters segments that lie below the mean (i.e. recently 
#' cultivated/harvested). If the segment length is greater than the 1st element in \emph{min.length} the 
#' segment is relabeled as "1 (i.e. "crop growth/maturity". This process is repeated for segments above 
#' the mean (i.e. crop growth/maturity). If the length of a segment is greater than the second element in 
#' \emph{min.length} it is labeled as "recently cultivated/harvested".}
#' @examples {
#' 
#' x <- c(293, 770, 1166, 1166, 1562, 2357, 3234, 
#' 5806, 5806, 5678, 5678, 5546, 5536, 5536, 5536, 
#' 5325, 5200, 4726, 3550, 2868, 2365, 2365, 2365)
#' 
#' n <- countCropCycles(x)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

countCropCycles <- function(x, min.length=c(1,1)) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if(!is.numeric(x)) {stop('"x" is not numeric')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. classify crop cylces
#---------------------------------------------------------------------------------------------------------------------#
  
  x1 <- x-mean(x) # identify break-point
  x1[x1 > 0] <- 1 # points above break-point (crop growth/maturity)
  x1[x1 < 0] <- 0 # points below break-point (recently cultivated/harvested)
  s = rle(as.numeric(x1)) # identify growth cycles
  
  # assign segment ID's (round I)
  s.id <- vector("numeric", length(x))
  for (i in 1:length(s$lengths)) {s.id[(sum(s$lengths[0:(i-1)])+1):sum(s$length[1:i])] <- i}
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. filter crop cylces
#---------------------------------------------------------------------------------------------------------------------#
  
  uv <- unique(s.id[x1==0]) # segments unique ID's (recently cultivated/harvested)
  vs <- sapply(uv, function(u){sum(s.id == u)})
  vi <- uv[which(vs < min.length[1])]
  for (u in 1:length(vi)) {x1[which(s.id == vi[u])] <- 1}
  
  uv <- unique(s.id[which(x1 == 1)]) # segments unique ID's (crop maturity)
  vs <- sapply(uv, function(u){sum(s.id == u)})
  vi <- uv[which(vs < min.length[2])]
  for (u in 1:length(vi)) {x1[which(s.id == vi[u])] <- 0}
  
  # assign segment ID's (round II)
  s.id <- vector("numeric", length(x))
  for (i in 1:length(s$lengths)) {s.id[(sum(s$lengths[0:(i-1)])+1):sum(s$length[1:i])] <- i}
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. count crop cylces
#---------------------------------------------------------------------------------------------------------------------#
  
  return(length(unique(s.id[x1==1])))
  
}

