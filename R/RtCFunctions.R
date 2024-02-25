## RESISTANCE TO CHANGE FUNCTIONS ##

# consyretr ---------------------------------------------------------------

#' Calculation consequence and feedback
#'
#' @description
#' @param wimp
#' @param std
#'
#' @return
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
