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
#'   indices based on the maximum total degree of the constructs. Defaults to FALSE
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
#' ph_indices_wnorm_non_std <- ph_index(wimp, method = "wnorm", std = FALSE)
#'

ph_index <- function(wimp, method = "weight", std = FALSE){

  # Connectivity of constructs
  c.io <- degree_index(wimp, method = method)

  # Rearrange In - Out columns
  c.io <- c.io[, c(2,1,3)]

  # Extract In and Out columns
  in.out <- c.io[, 1:2]

  # Linear Transformation matrix
  coef <- 0.5 * sqrt(2)
  coef.matrix <- matrix(c(coef, -coef, coef, coef), nrow = 2)

  # Calculate P - H matrix
  ph.mat <- in.out %*% t(coef.matrix)
  colnames(ph.mat) <- c("p", "h")

  # Standardization
  if (std){
    max.total.deg <- max(c.io[, 3])
    ph.mat <- ph.mat / max.total.deg
  }

  return(ph.mat)
}

# Mahalanobis Indices ------------------------------------------------------------

#' Calculate Mahalanobis Indices for Constructs
#'
#' This function calculates the Mahalanobis distance for each construct
#' in a given PH matrix obtained from a `wimp` object. It also determines
#' whether each construct is considered "central" based on a chi-square
#' cutoff, which is calculated using a specified significance level.
#'
#' @param wimp A `wimp` object containing the data from which the PH matrix
#'   is derived.
#' @param method A character string specifying the method used for calculating
#'   the PH matrix. Default is `"weight"`.
#' @param std A logical value indicating whether the data should be
#'   standardized before calculating the Mahalanobis distance. Default is `FALSE`.
#'
#' @return A matrix that includes the original PH matrix with two additional
#'   columns: one for the Mahalanobis distance of each construct and another
#'   indicating whether the construct is considered "central" based on the
#'   chi-square cutoff.
#'
#' @export
#'
#' @examples
#' result <- mahalanobis_index(wimp, method = "weight", std = FALSE)
#' head(result)
#'

mahalanobis_index <- function(wimp, method = "weight", std = FALSE){

  # Obtain PH matrix associated to the wimp object
  ph.mat <- ph_index(wimp = wimp, method = method, std = std)

  # Mahalanobis distance--------

  # Covariance matrix
  cov.matrix <- cov(ph.mat)
  # Means vector
  means.vect <- colMeans(ph.mat)
  # Calculate the Mahalanobis distance for each observation in the dataframe
  ph.mat.df <- as.data.frame(ph.mat)
  m.dist <- mahalanobis(ph.mat.df, center = means.vect, cov = cov.matrix)

  # Add column to the results matrix
  phm.mat <- cbind(ph.mat, m.dist)

  # Determine the chi-square cutoff for the given significance level--------------
  sign.level <- 0.2
  # Degrees of freedom
  df <- ncol(ph.mat)
  # Chi-Square Cutoff
  chi.square.cutoff <- qchisq(1 - sign.level, df)

  # Annotate observations as "central" contructs based on the chi-square cutoff
  central <- ifelse(phm.mat[,"m.dist"] > chi.square.cutoff, TRUE, FALSE)

  # Add column to the results matrix
  phmc.mat <- cbind(phm.mat, central)
  return(phmc.mat)
}
