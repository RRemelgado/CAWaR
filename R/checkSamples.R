#' @title checkSamples
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description checks if a shapefile with ground truth data contains all the necessary fields in the correct format.
#' @param x An object or a list of class \emph{sp} containing a \emph{data.frame} (e.g. \emph{SpatialPolygonsDataFrame}).
#' @return A \emph{data.frame} with the consistency checks for each element in \emph{x}.
#' @importFrom raster shapefile
#' @importFrom sp coordinates
#' @details {Checks if a shapefile - or a list of - contains necessary columns and if these have the right format. It searches for:
#' \itemize{
#'  \item{\emph{sampler} - Character vector with name of responsible person.}
#'  \item{\emph{date} - Date vector with the date on which each sample was collected (formated as "yyyy-mm-dd").}
#'  \item{\emph{label} - Character vector sample label (e.g. land cover class).}}}
#' @seealso \code{\link[fieldRS]{labelCheck}}
#' @examples {
#' 
#' require(fieldRS)
#' 
#' # Example ground-truth data
#' data(fieldData)
#' 
#' # check shapefile content
#' cs <- checkSamples(fieldData)
#' head(cs)
#' 
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

checkSamples <- function(x) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if(!is.list(x)) {x <- list(x)}
  c <- sapply(x, function(s) {return(ifelse(class(try(s@data, silent=TRUE))[1]=="try-error", FALSE, TRUE))})
  if (sum(c)!=length(c)) {
    warning('one or more elements is "x" are a valid spatial object (check the function output for a failure index)')
    return(c)}
  
#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. find errors
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  control <- lapply(x, function(s) {
    
    col.names <- tolower(colnames(s@data)) # list column names
    
    # check date field (when was the sample collected?)
    if (!'date' %in% col.names) {date.test <- 'ERROR (missing "date" column)'} else {
      date.test <- try(as.Date(s@data$date), silent=TRUE)
      if (class(date.test)[1] == 'try-error') {date.test <- 'ERROR (Not a Date object)'} else {
        if (sum(is.na(date.test))==0) {date.test <- 'PASSED'} else {date.test <- 'NOTE (found missing values)'}}}
    
    # check sampler field (who was responsible?)
    if (!'sampler' %in% col.names) {sampler.test <- 'ERROR (missing "sampler" column)'} else {
      if (!is.character(s@data$sampler)) {sampler.test <- 'ERROR ("sampler" is not a character vector)'} else {
        if (sum(is.na(s@data$sampler))==0) {sampler.test <- 'PASSED'} else {sampler.test <- 'NOTE (found missing values)'}}}
    
    # check label field (what is the sample from?)
    if (!'label' %in% col.names) {label.test <- 'ERROR (missing "label" column)'} else {
      if (!is.character(s@data$label)) {label.test <- 'ERROR ("label" is not a character vector)'} else {
        if (sum(is.na(s@data$label))==0) {label.test <- 'PASSED'} else {sampler.test <- 'NOTE (found missing values)'}}}
    
    # compile test results
    control <- data.frame(sampler=sampler.test, date=date.test, label=label.test, stringsAsFactors=FALSE)
    return(control)

  })

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. provide final report
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  control <- do.call(rbind, control)
  return(control)

}
