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
#' density_index(example_wimp)

density_index <- function(wimp) {

  # Get weights matrix directly from wimp
  if (!is.null(wimp$global$weight_matrix)) {
    wmat <- wimp$global$weight_matrix
  } else {
    # Fallback: reconstruct from edges
    edges <- wimp$edges
    n <- nrow(wimp$vertices)
    wmat <- matrix(0, nrow = n, ncol = n)
    if (!is.null(edges) && nrow(edges) > 0) {
      for (r in seq_len(nrow(edges))) {
        i <- as.integer(edges[r, "from"])
        j <- as.integer(edges[r, "to"])
        wmat[i, j] <- as.numeric(edges[r, "weight"])
      }
    }
  }
  n <- ncol(wmat)
  
  # Calculate density as E / (n * (n-1))
  # Count edges using 'simple' method in degree_index
  result <- sum(degree_index(wimp, method = "simple")[, 1]) / (n * (n - 1))
  
  result
}
