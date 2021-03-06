#' @title phenoCropVal
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Spatially explicit and phenology driven validation scheme for cropland mapping.
#' @param x A \emph{matrix} or \emph{data.frame}.
#' @param y A \emph{character} vector.
#' @param z A \emph{character} vector.
#' @return A \emph{list} containing a set of reference profiles for each unique class in \emph{y}.
#' @importFrom stats cor
#' @importFrom raster which.max
#' @importFrom ggplot2 ggplot aes_string geom_bar
#' @details {For each unique class in \emph{y}, the function iterates through each unique element in \emph{z} 
#' and keeps it for validation. Then, it calls \code{\link{analyseTS}} to derive reference profiles for each 
#' unique class in \emph{y} and uses them to classify the validation samples using \code{\link{phenoCropClass}}. 
#' The final output consists of:
#' \itemize{
#'  \item{\emph{sample.validation} - A \emph{logical} vector with the same length of \emph{x} where TRUE means it was correctly classified.}
#'  \item{\emph{predicted.class} - A \emph{character} vector with the predicted classes for each sample.}
#'  \item{\emph{sample.count} - A \emph{numeric} vector with the number of non-NA used for validation per sample.}
#'  \item{\emph{sample.r2} - A \emph{numeric} vector with the r2 value between the target sample and the selected class profile.}
#'  \item{\emph{class.accuracy} - A \emph{data.frame} with sample count per class, precision, recall and F1-scores per unique class in \emph{y}.}}}
#' @seealso \code{\link{extractTS}} \code{\link{phenoCropClass}}
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
#' 
#' # read reference profiles
#' data(referenceProfiles)
#' 
#' # read time series
#' data(fieldDataTS)
#' fieldDataTS <- as.data.frame(fieldDataTS$weighted.mean)
#' 
#' # read info. on sample spatial grouping
#' data(fieldDataCluster)
#' 
#' # derive validation results
#' cropVal <- phenoCropVal(fieldDataTS, fieldData$crop, fieldDataCluster$region.id)
#' 
#' # plot accuracy results
#' cropVal$accuracy.plot
#' 
#' # plot correctly classified polygons in red
#' plot(fieldData)
#' plot(fieldData[cropVal$sample.validation,], col="red", add=TRUE)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

phenoCropVal <- function(x, y, z) {
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. Check variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  if (class(x)[1] != "data.frame") {stop('"x" is not a data.frame')}
  
  if (missing(y)) {stop('"y" is missing')}
  if (!is.character(y)) {stop('"y" is not a character vector')}
  if (nrow(x) != length(y)) {stop('"x" and "y" have different lenghts')}
  unique.y <- unique(y)
  
  if (missing("z")) {
    warning('"z" is missing (each entry in "x" will be validated separately)')
    z <- 1:length(y)
  } else {if (length(y) != length(z)) {stop('"y" and "z" have different lenghts')}}
  
  i0 <- 1:nrow(x)
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. Valiate samples and estimate class accuracy
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  sample.val <- vector("logical", nrow(x)) # sample-wise validation
  sample.class <- vector("character", nrow(x))
  class.sum <- sample.sum <- sample.true <- vector("numeric", length(unique.y))
  sample.count <- vector("numeric", nrow(x))
  sample.r2 <- vector("numeric", nrow(x))
  
  for (c in 1:length(unique.y)) {
    
    unique.z <- unique(z[which(y == unique.y[c])]) # indices of target samples
    
    if (length(unique.z) > 1) {
      
      for (k in 1:length(unique.z)) {
        
        # determine validation/training indices
        vi <- which(z == unique.z[k] & y == unique.y[c])
        ti <- i0[!i0 %in% vi]
        
        # assign class
        cropClass <- do.call(rbind, lapply(1:length(vi), function(j) {
          pc <- phenoCropClass(as.numeric(x[vi[j],]), x[ti,])
          i <- which.max(pc$r2)[1]
          pc <- pc[i,]
          pc$class <- y[ti[i]]
          return(pc)}))
        
        # 
        sample.class[vi] <- cropClass$class
        sample.count[vi] <- cropClass$count
        sample.r2[vi] <- cropClass$r2
        
        # validate results
        sample.val[vi] <- y[vi] == sample.class[vi]
        
        sample.true[c] <- sample.true[c] + sum(sample.val[vi] & y[vi]==unique.y[c]) # count class occurrences
        class.sum <- class.sum + sapply(unique.y, function(u) {sum(sample.class[vi] == u)}) # count of class occurrence
        sample.sum[c] <- sample.sum[c] + length(vi) # number of validation samples
        
      }
      
    }
    
  }
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 3. Derive accuracies
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  p <- sample.true / sample.sum
  r <- sample.true / class.sum
  class.acc <- data.frame(class=unique.y, count=sample.sum, precision=p, recall=r, f1=(2 * ((p * r) / (p + r))), stringsAsFactors=FALSE)
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 4. Build accuracy plot
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  p <- ggplot(class.acc, aes_string(x="class", y="f1")) + geom_bar(stat="identity")
  
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 5. Derive output
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(sample.validation=sample.val, predicted.class=sample.class, sample.ts.count=sample.count, 
              sample.r2=sample.r2, class.accuracy=class.acc, accuracy.plot=p))
  
}
