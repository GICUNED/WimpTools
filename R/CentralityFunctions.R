## CENTRALITY FUNCTIONS ##

# PH Indices ------------------------------------------------------------

#' Title
#'
#' @description To obtain the presence (P, frequency of occurrence) and hierarchy
#' (H, influence on others) of constructs associated with an implication grid.
#'
#' @param wimp
#'
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#'

index_ph <- function(wimp, method = "simple"){

  # Connectivity of constructs
  c.io <- degree_index(wimp, method = method)

  # Linear transformation matrix
  coef <- 0.5 * sqrt(2)
  mat.proy <- matrix(c(coef, coef, -coef, coef), nrow = 2, byrow = TRUE)

  # Extract In and Out columns
  in.out <- c.io[, 1:2]

  # Calculate p and h for each row
  p <- rowSums(in.out * c(coef, coef))
  h <- rowSums(in.out * c(-coef, coef))

  # Standardization
  max.total.deg <- max(c.io[, 3])
  p <- p / max.total.deg
  h <- h / max.total.deg

  # Create the new matrix with p and h dimensions
  ph.mat <- matrix(c(p, h), ncol = 2, byrow = FALSE,
                   dimnames = list(rownames(c.io), c("p", "h")))

  return(ph.mat)
}

