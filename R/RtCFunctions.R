## RESISTANCE TO CHANGE FUNCTIONS ##

# consyretr ---------------------------------------------------------------

#' Calculation of consequence and feedback for each construct
#'
#' @description This function calculates the number of consequences and feebacks
#' for each construct to analyze resistance to change
#'
#' @param wimp
#' @param std
#'
#' @return Returns a dataframe with number of consequences and feedbacks
#'
#' @import
#'
#' @export
#'
#' @examples
#' cyr_index (example.wimp)
#'
cyr_index <- function(wimp, std=FALSE) {
  wimp <- .align.wimp(wimp)
  wmatrix<- wimp $scores $weights
  consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
  feedback <- wmatrix * t (wmatrix)
  feedback <- apply (feedback, MARGIN = 1, FUN = sum)
  result <- data.frame (feedback , consequence)
  return(result)
  }
