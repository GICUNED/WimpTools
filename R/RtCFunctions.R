## RESISTANCE TO CHANGE FUNCTIONS ##

# cq.fb ---------------------------------------------------------------

#' Calculation of consequence and feedback effects for each construct of the
#' system
#'
#' @description This function calculates the numbers of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp
#' @param std
#'
#' @return Returns a dataframe with numbers of consequences and feedbacks
#'
#' @import
#'
#' @export
#'
#' @examples
#' cq.fb_index (example.wimp)
#'
cq.fb_index <- function(wimp, std=FALSE) {
  wimp <- .align.wimp(wimp)
  wmatrix<- wimp $scores $weights
  consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
  feedback <- wmatrix * t (wmatrix)
  feedback <- apply (feedback, MARGIN = 1, FUN = sum)
  result <- data.frame (feedback , consequence)
  return(result)
  }
