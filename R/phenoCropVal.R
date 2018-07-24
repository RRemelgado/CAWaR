#' @title phenoCropVal
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Spatialy explicit and phenology driven validation scheme for cropland mapping.
#' @param x A \emph{matrix} or \emph{data.frame}.
#' @param y A \emph{character} vector.
#' @param z A \emph{character} vector.
#' @return A \emph{list} containing a set of reference profiles for each unique class in \emph{y}.
#' @importFrom stats cor
#' @importFrom raster which.max
#' @importFrom ggplot2 ggplot aes_string geom_bar
#' @details {For each unique class in \emph{y}, the function iterates through each unique element in \emph{z} 
#' and keeps it for validation. Then, it calls \code{\link{analyzeTS}} to derive reference profiles for each 
#' unique class in \emph{y} and uses them to classify the validation samples using \code{\link{phenoCropClass}}. The 
#' final output consists of:
#' \itemize{
#'  \item{\emph{sample.validation} - A \emph{logical} vector with the same length of \emph{x} where TRUE means it was correctly classied.}
#'  \item{\emph{class.accuracy} - A \emph{data.frame} with an F1-score for each unique class in \emph{y}.}}}
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
#' # derive time series
#' ev <- as.data.frame(extractTS(fieldData, extend(r, 60))$weighted.mean)
#' 
#' # derive validation results
#' cropVal <- phenoCropVal(ev, fieldData$crop)
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
    warning('"z" is missing (each entry in "x" will be validated separately')
    z <- 1:length(y)
  } else {if (length(y) != length(z)) {stop('"y" and "z" have different lenghts')}}
  
  i0 <- 1:nrow(x)
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Valiate samples and estimate class accuracy
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  sample.val <- vector("logical", nrow(x)) # sample-wise validation
  class.acc <- data.frame(class=unique.y, f1=0, stringsAsFactors=FALSE) # class-wise validation
  
  for (c in 1:length(unique.y)) {
    
    unique.z <- unique(z[which(y == unique.y[c])]) # indices of target samples
    
    pp <- 0
    cp <- 0
    ss <- 0
    
    for (k in 1:length(unique.z)) {
      
      # determine validation/training indices
      vi <- which(z == unique.z[k] & y == unique.y[c])
      ti <- i0[!i0 %in% vi]
      
      # build reference profiles
      tmp <- analyzeTS(x[ti,], y[ti])
      reference.ts <- do.call(rbind, lapply(tmp$y.statistics, function(j) {j$median}))
      reference.class <- tmp$labels
      
      # validate results
      sample.val[vi] <- y[vi] == reference.class[sapply(1:length(vi), function(j) {which.max(phenoCropClass(as.numeric(x[vi[j],]), reference.ts, 1)$r2)})]
      
      pp <- pp + sum(sample.val[vi]) # count class occurrences
      cp <- cp + sum(sample.val[vi] & y[vi]==unique.y[c]) # count correct occurrences
      ss <- ss + length(vi) # number of validation samples
      
    }
    
    # estimate F1-score
    p <- cp / pp
    r <- cp / ss
    class.acc$f1[c] <- 2 * ((p * r) / (p + r))
    
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Build accuracy plot
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  p <- ggplot(class.acc, aes_string(x="class", y="f1")) + geom_bar(stat="identity")
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Derive output
#-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  return(list(sample.validation=sample.val, class.accuracy=class.acc, accuracy.plot=p))
  
}
