## CENTRALITY FUNCTIONS ##

# PH Indices ------------------------------------------------------------

#' Calculate Presence and Hierarchy Indices
#'
#' This function computes the presence (P, frequency of occurrence) and
#' hierarchy (H, influence on others) indices for constructs within an implication grid.
#' It can standardize these indices based on the maximum degree if required.
#'
#' @param wimp An object of class 'wimp', which contains an implication grid
#'   and associated constructs.
#' @param method A character string specifying the method used to calculate
#'   the degree indices. Default is "weight". Other methods may be available
#'   depending on the implementation of 'degree_index' function.
#' @param std A logical value indicating whether to standardize the P and H
#'   indices based on the maximum total degree of the constructs. Defaults to TRUE.
#'
#' @return A 2-column matrix with the presence index 'p' and the hierarchy index 'h'
#'   for each construct. If 'std' is TRUE, these values are standardized.
#'
#' @export
#'
#' @examples
#' # Assuming 'wimp' is an object with an implication grid
#' ph_indices <- ph_index(wimp)
#' ph_indices_std <- ph_index(wimp, std = TRUE)
#' ph_indices_non_std <- ph_index(wimp, std = FALSE)
#'

ph_index <- function(wimp, method = "weight", std = TRUE){

  # Connectivity of constructs
  c.io <- degree_index(wimp, method = method)

  # Rearrange In - Out columns
  c.io <- c[, c(2,1,3)]
  # Extract In and Out columns
  in.out <- c.io[, 1:2]

  # Calculate p and h for each row
  coef <- 0.5 * sqrt(2)
  p <- rowSums(c(coef, coef) * in.out)
  h <- rowSums(c(-coef, coef) * in.out)

  # Standardization
  if (std){
    max.total.deg <- max(c.io[, 3])
    p <- p / max.total.deg
    h <- h / max.total.deg
  }

  # Create the new matrix with p and h dimensions
  ph.mat <- matrix(c(p, h), ncol = 2, byrow = FALSE,
                   dimnames = list(rownames(c.io), c("p", "h")))

  return(ph.mat)
}

