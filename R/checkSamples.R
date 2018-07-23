#' @title checkSamples
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description checks if a shapefile with ground truth data contains all the necessary fields in the correct format.
#' @param samples Vector of class \emph{character}.
#' @return A \emph{data.frame} with the consistency checks for each element in \emph{samples}.
#' @importFrom raster shapefile
#' @details {Checks if a shapefile - or a set of - contains necessary columns and if these have the right format. It searches for:
#' \itemize{
#'  \item{\emph{sampler} - Character vector with name of responsible person.}
#'  \item{\emph{date} - Date vector with the date on which each sample was collected (formated as "yyyy-mm-dd").}
#'  \item{\emph{label} - Character vector sample label (e.g. land cover class).}}}
#' @seealso \code{\link[fieldRS]{labelCheck}}
#' @examples {}
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

checkSamples <- function(samples=samples) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (!is.character(samples)) {stop('"samples" is not a character vector')}

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. find errors
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  control <- lapply(samples, function(s) {

    if (!file.exists(s)) {control <- 'File not found'} else {

      file <- try(shapefile(s), silent=TRUE) # check if file can be read

      if (class(file)[1] != 'SpatialPolygonsDataFrame') {control <- 'file is not a "SpatialPolygonssDataFrame"'} else {

        col.names <- tolower(colnames(file@data)) # list column names

        # check date field (when was the sample collected?)
        if (!'date' %in% col.names) {date.test <- 'ERROR (missing "date" column)'} else {
          date.test <- try(as.Date(file@data$date), silent=TRUE)
          if (class(date.test)[1] == 'try-error') {date.test <- 'ERROR (Not a Date object)'} else {
            if (sum(is.na(date.test))==0) {date.test <- 'PASSED'} else {date.test <- 'NOTE (found missing values)'}}}

        # check sampler field (who was responsible?)
        if (!'sampler' %in% col.names) {sampler.test <- 'ERROR (missing "sampler" column)'} else {
          if (!is.character(file@data$sampler)) {sampler.test <- 'ERROR ("sampler" is not a character vector)'} else {
            if (sum(is.na(file@data$sampler))==0) {sampler.test <- 'PASSED'} else {sampler.test <- 'NOTE (found missing values)'}}}

        # check label field (what is the sample from?)
        if (!'label' %in% col.names) {label.test <- 'ERROR (missing "label" column)'} else {
          if (!is.character(file@data$label)) {label.test <- 'ERROR ("label" is not a character vector)'} else {
            if (sum(is.na(file@data$label))==0) {label.test <- 'PASSED'} else {sampler.test <- 'NOTE (found missing values)'}}}

        # compile test results
        control <- data.frame(sampler=sampler.test, date=date.test, label=label.test, stringsAsFactors=FALSE)
        return(control)

      }}})

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. provide final report
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  control <- do.call(rbind, control)
  row.names(control) <- basename(samples)
  return(control)

}
