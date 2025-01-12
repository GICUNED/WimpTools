## WIMPGRID INDICES FUNCTIONS ##

#' Digraph Density Index -- density_index()
#'
#' @description Function used to calculate the density of edges of the
#' calculated digraph of the impgrid
#'
#' @param wimp  Subject's Weigthed ImpGrid. It must be a "wimp" S3 object
#' imported by the \code{\link{importwimp}} function.
#'
#' @return Returns a value from 0 to 1 representing the ratio of the number of
#' edges in the graph over the maximum number of possible edges.
#'
#' @export
#'
#' @examples
#'
#' density_index(example.wimp)

density_index <- function(wimp) {

  wmat <- wimp$scores[[3]]
  n <- ncol(wmat)

  result <- sum(degree_index(wimp)[,1]) / (n * (n - 1))

  return(result)
}
