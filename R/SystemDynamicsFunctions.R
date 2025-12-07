## System Dynamics Functions ##

# Scenario Matrix -----------------------------------------------------------

#' Scenario Matrix -- scenariomatrix()
#'
#' @description Generate a scenario matrix from a wimp object and action
#'   vector. Iteratively simulates outcomes using the selected inference
#'   method to support scenario analysis.
#'
#' @param wimp A "wimp" S3 object.
#' @param infer Inference type: "self dynamics" or "impact dynamics".
#'   Defaults to "self dynamics".
#' @param thr Thresholding method: "saturation", "tanh", or "none".
#'   Defaults to "saturation".
#' @param act_vector Numeric vector of construct activations applied to the
#'   grid. Length must match number of constructs unless using impact dynamics.
#' @param max_iter Maximum number of iterations. Defaults to 5.
#' @param e Convergence tolerance. Smaller values are stricter. Default 0.0001.
#' @param stop_iter Consecutive iterations below tolerance required to stop.
#'   Defaults to 2.
#' @param exclude_dilemmatics If TRUE, exclude dilemmatic constructs in
#'   calculations. Defaults to TRUE.
#'
#' @return A "scn" S3 class object.
#' @export
#'
#' @examples
#'
#' # Activation vector (example)
#' act_vector <- c(0.1, -0.2, 0.3, -0.1, 0)
#'
#' # Generate a scenario matrix
#' scenariomatrix(wimp = example_wimp, act_vector = act_vector)
#'

scenariomatrix <- function(wimp, act_vector = NA, infer = "self dynamics",
                           thr = "saturation", max_iter = 5, e = 0.0001,
                           stop_iter = 2, exclude_dilemmatics = FALSE) {

  if (!inherits(wimp, "wimp")) {
    stop("The Weighted Implication Grid must be class wimp.")
  }

  # Validate weight matrix and activation vector dimensions
  if (is.null(wimp$global$weight_matrix)) {
    stop("wimp object must contain weight_matrix in wimp$global")
  }
  act_vector_input <- act_vector
  n_constructs <- nrow(wimp$global$weight_matrix)
  if (length(act_vector) != n_constructs && infer != "impact dynamics") {
    stop("Length of act_vector (", length(act_vector),
         ") must match number of constructs (", n_constructs, ").")
  }

  # Align constructs so ideal values are positive
  wimp <- .align_wimp(wimp, exclude_dilemmatics = exclude_dilemmatics)

  # Extract ideal values for swap vector
  ideal <- wimp$vertices$ideal
  swap_vector <- ideal / abs(ideal)
  swap_vector[is.nan(swap_vector)] <- 1

  if (infer == "self dynamics") {
    act_vector <- act_vector * swap_vector
  }

  # Initialize scenario matrix from self values
  n_constructs <- nrow(wimp$vertices)
  scene_matrix <- t(matrix(wimp$vertices$self))
  trans_matrix <- t(wimp$global$weight_matrix)
  next_matrix <- trans_matrix

  n <- 1
  i <- 0

  # Iterate to convergence or max iterations
  while (n <= max_iter && i <= stop_iter) {

    if (infer == "self dynamics") {
      # Calculate next iteration with activation effects
      next_iter <- scene_matrix[n, ] + t(act_vector)
      next_iter <- mapply(.thr, next_iter, thr)

      # Calculate delta and update scenario matrix
      delta_iter <- next_iter - scene_matrix[n, ]
      scene_matrix <- rbind(scene_matrix, next_iter)

      # Update activation vector for next iteration
      act_vector <- trans_matrix %*% delta_iter
    }

    if (infer == "impact dynamics") {
      # Initialize impact scenario matrix on first iteration
      if (n == 1) {
        scene_matrix <- t(rep(0, n_constructs))
      }

      # Set dilemmatic constructs to zero if requested
      if (!exclude_dilemmatics) {
        n_matrix <- next_matrix
        n_matrix[.which_dilemmatics(wimp), ] <- 0
      } else {
        n_matrix <- next_matrix
      }

      # Calculate impact totals
      sum_columns <- t(n_matrix) %*% rep(1, nrow(trans_matrix))
      next_iter <- t(sum_columns)
      scene_matrix <- rbind(scene_matrix, next_iter)
      next_matrix <- trans_matrix %*% next_matrix
    }

    # Check convergence criterion
    e_iter <- mean(abs(next_iter - scene_matrix[n, ]))

    if (e_iter < e) {
      i <- i + 1
    } else {
      i <- 0
    }
    n <- n + 1
  }

  # Name rows of scenario matrix with iteration numbers
  rownames(scene_matrix) <- paste("iter", 0:(n - 1))
  colnames(scene_matrix) <- paste(wimp$vertices$left_pole, " - ",
                                  wimp$vertices$right_pole, sep = "")

  # Check convergence status
  if (n < max_iter) {
    convergence <- n - (stop_iter + 1)
  } else {
    convergence <- NA
  }

  # Build scenario list for S3 scn class
  scene_list <- list()
  scene_list$values <- scene_matrix
  scene_list$convergence <- convergence

  # Build constructs data structure for scn object
  scene_list$constructs <- list(
    left_pole = wimp$vertices$left_pole,
    right_pole = wimp$vertices$right_pole,
    constructs = paste(wimp$vertices$left_pole, " - ",
                       wimp$vertices$right_pole, sep = "")
  )
  scene_list$self <- list(
    self = wimp$vertices$self,
    ideal = wimp$vertices$ideal
  )
  scene_list$weights <- wimp$global$weight_matrix
  scene_list$method <- list(
    infer = infer,
    threshold = thr
  )
  scene_list$params <- list(
    infer = infer,
    threshold = thr,
    max_iter = max_iter,
    e = e,
    stop_iter = stop_iter,
    exclude_dilemmatics = exclude_dilemmatics,
    act_vector = act_vector_input
  )

  class(scene_list) <- c("scn", "list")

  return(scene_list)
}

# PCSD -----------------------------------------------------------------

#' Personal Constructs System Dynamics plot -- pcsd()
#'
#' @description Interactive line plot of personal constructs system dynamics.
#'   Displays \code{\link{scenariomatrix}} values across iterations.
#'
#' @param scn A "scn" object returned by \code{scenariomatrix()}.
#' @param vline Optional iteration index to highlight with a vertical line.
#'   Defaults to NA (no line).
#'
#' @return Interactive plot created with plotly.
#'
#' @import plotly
#'
#' @export
#'
#' @examples
#'
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example_wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act_vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Plot the dynamics
#' pcsd(example_scn)
#'


pcsd <- function(scn, vline = NA) {
  # Extract constructs and configuration
  poles <- scn$constructs$constructs
  n_constructs <- length(poles)
  infer <- scn$method$infer
  iter <- nrow(scn$values)

  # Prepare self matrix and result values
  self_vector <- scn$self$self
  self_matrix <- matrix(self_vector, ncol = length(self_vector),
                        nrow = iter, byrow = TRUE)
  res <- scn$values

  # Create dataframe based on inference type
  x <- 0:(iter - 1)
  if (infer == "self dynamics") {
    df <- data.frame(x, res - self_matrix)
  }
  if (infer == "impact dynamics") {
    df <- data.frame(x, res / n_constructs)
  }

  # Rename columns for easier reshaping
  colnames(df) <- c("x", poles)

  # Convert to long format for Plotly
  df_long <- tidyr::pivot_longer(
    df, cols = -x, names_to = "construct", values_to = "value"
  )

  # Define dynamic Y axis range
  maxv <- max(abs(df_long$value)) * 1.05

  # Plot using automatic color and symbol mapping
  fig <- plot_ly(
    data = df_long,
    x = ~x,
    y = ~value,
    color = ~construct,
    symbol = ~construct,
    symbols = c("circle", "square", "diamond", "cross", "x",
                "triangle-up", "triangle-down", "triangle-left",
                "triangle-right", "star"),
    type = "scatter",
    mode = "lines+markers",
    marker = list(size = 10),
    line = list(shape = "spline")
  ) %>%
    layout(
      xaxis = list(
        title = list(
          text = "ITERATIONS",
          font = list(size = 20)
        ),
        tickfont = list(size = 20)
      ),
      yaxis = list(
        title = list(
          text = "SELF DIFFERENTIAL",
          font = list(size = 20)
        ),
        tickfont = list(size = 20),
        range = c(-maxv, maxv)
      ),
      legend = list(
        title = list(
          text = "<b>PERSONAL CONSTRUCTS</b>",
          font = list(size = 15)
        ),
        font = list(size = 20)
      )
    )

  # Add vertical reference line if specified
  if (!is.na(vline)) {
    fig <- fig %>% add_lines(
      x = vline,
      y = c(-maxv, maxv),
      line = list(color = "grey", dash = "dot"),
      inherit = FALSE,
      showlegend = FALSE
    )
  }

  return(fig)
}


# AUC Index ---------------------------------------------------------------

#' PCSD AUC Index -- auc_index()
#'
#' @description This function calculates the area under the PCSD curve for each
#' construct.
#'
#' @param scn A "scn" object returned by \code{scenariomatrix()}.
#' @param std If \code{TRUE}, standardizes the AUC results. Defaults to TRUE.
#'
#' @return Returns a vector with the AUC index of each construct.
#'
#' @import MESS
#'
#' @export
#'
#' @examples
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example_wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act_vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate AUC index for each construct
#' auc_index(example_scn, std = TRUE)
#'

auc_index <- function(scn, std = TRUE) {

  poles <- scn$constructs$constructs
  iter <- nrow(scn$values)
  self_vector <- scn$self$self
  self_matrix <- matrix(
    self_vector, ncol = length(self_vector), nrow = iter, byrow = TRUE
  )
  std_coef <- ifelse(std, iter, 1)

  res <- scn$values

  if (scn$method$infer == "self dynamics") {
    res <- res - self_matrix
  }

  auc_matrix <- matrix(ncol = length(poles), nrow = 1)

  for (n in seq_along(poles)) {
    auc_matrix[, n] <- MESS::auc(
      c(0:(iter - 1)), res[, n], type = "spline"
    ) / std_coef
  }

  result <- t(auc_matrix)

  rownames(result) <- poles
  colnames(result) <- "AUC"

  return(result)
}

# PCSD Stability Index ----------------------------------------------------

#' PCSD Stability Index -- stability_index()
#'
#' @description Compute the standard deviation of each construct across PCSD
#'   iterations.
#'
#' @param scn A "scn" object returned by \code{scenariomatrix()}.
#'
#' @return A vector with the standard deviation of each construct.
#'
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom tidyr pivot_longer
#'
#' @export
#'
#' @examples
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example_wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act_vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate stability index for each construct
#' stability_index(example_scn)
#'

stability_index <- function(scn) {

  poles <- scn$constructs$constructs
  res <- scn$values

  # Standard deviation for each construct across iterations
  result <- apply(res, 2, sd)

  result <- matrix(result)
  rownames(result) <- poles
  colnames(result) <- "Standard Deviation"

  result
}

# PCSD Summary ------------------------------------------------------------

#' PCSD summary -- pcsd_summary()
#'
#' @description Summarize PCSD results: initial, final, average, and change per
#'   construct.
#'
#' @param scn A "scn" object returned by \code{scenariomatrix()}.
#'
#' @return A matrix with PCSD summary metrics per construct.
#'
#'
#' @export
#'
#' @examples
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example_wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act_vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate PCSD summary for each construct
#' pcsd_summary(example_scn)
#'

pcsd_summary <- function(scn) {

  poles <- scn$constructs$constructs
  iter <- nrow(scn$values)

  res <- scn$values
  average <- colMeans(res)

  # First and last iterations (transpose to constructs x metrics)
  result <- res[c(1, iter), ]
  result <- t(result)
  result <- cbind(result, average, result[, 2] - result[, 1])

  rownames(result) <- poles
  colnames(result) <- c("Initial", "Final", "Average", "Difference")

  result
}
