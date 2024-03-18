## RESISTANCE TO CHANGE FUNCTIONS ##

# cq.fb ---------------------------------------------------------------

#' Calculation of consequence and feedback effects for each construct of the
#' system
#'
#' @description This function calculates the intensity of the effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp
#' @param std
#'
#' @return Returns a dataframe with the intensity of the effects of consequences and feedbacks
#'
#'
#'
#' @export
#'
#' @examples
#' cq.fb_index (example.wimp)
#'
cq.fb_index <- function(wimp, std=FALSE) {
  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights
  # Extract punctuations of consequences as measured when comparing to the hypothetical self
  consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
  # Calculate punctuations of feedback effects by multiplying the weight matrix by its own transposed version
  feedback <- wmatrix * t (wmatrix)
  # Sum feedback punctuations to get a feedback index for each construct polarity
  feedback <- apply (feedback, MARGIN = 1, FUN = sum)
  # Define the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (feedback , consequence)
  # Setting names to the poles of each conastruct
  poles <- wimp[["constructs"]][["constructs"]]
  # Show row names (titles)
  rownames(result) <- poles
  return(result)
}

# rtc_plot ----------------------------------------------------------------

#' Visualization of the effects of Consequences and Feedback
#'
#' @description This function creates a scatter plot to show the effects of the consequences and feedbacks of each construct to analyze their resistance to change
#'
#' @param wimp
#' @param std
#'
#' @return returns a scatter plot showing the effects of consequences and feedback on each construct of the system
#' @export
#' @import plotly
#' @examples
#' rtc_plot (example.wimp)
rtc_plot <- function(wimp, std=FALSE) {
  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights
  # Extract punctuations of consequences as measured when comparing to the hypothetical self
  consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
  # Calculate punctuations of feedback effects by multiplying the weight matrix by its own transposed version
  feedback <- wmatrix * t (wmatrix)
  # Sum feedback punctuations to get a feedback index for each construct polarity
  feedback <- apply (feedback, MARGIN = 1, FUN = sum)
  # Define the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (feedback , consequence)
  # Setting names to the poles of each conastruct
  poles <- wimp[["constructs"]][["constructs"]]
  # Show row names (titles)
  rownames(result) <- poles

# Create a scatter plot to visualize the effects of cq and fb
  fig <- plot_ly(data = result, x = ~consequence, y = ~feedback, color = ~poles)
  return(fig)
}




