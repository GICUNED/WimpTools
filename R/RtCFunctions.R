## RESISTANCE TO CHANGE FUNCTIONS ##

# cq.fb ---------------------------------------------------------------

#' Calculation of consequence and feedback effects for each construct of the
#' system
#'
#' @description This function calculates the numbers of consequences and feedbacks
#' @description This function calculates the intensity of the effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp
#' @param std
#' @param std A logical value indicating wether to standarize the scores (TRUE) or (FALSE)
#'
#' @return Returns a dataframe with numbers of consequences and feedbacks
#' @return Returns a dataframe with the intensity of the effects of consequences and feedbacks
#'
#'
#'
#' @export
#'
#' @examples
#' cq.fb_index (example.wimp)
#'

cq.fb_index <- function(wimp, std) {

  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights

  if (std == FALSE) {
    # Extract scores of consequences as measured when comparing to the hypothetical self
    consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
    # Calculate scores of feedback effects by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)
    # Extract feedback punctuations to get a feedback index for each construct polarity
    feedback <- apply (feedback, MARGIN = 1, FUN = sum)

  } else {
    # Extract normalized scores (divided by the column's mean) of consequences as measured when comparing to the hypothetical self
    consequence <- apply (wmatrix, MARGIN = 1, FUN = mean)
    # Calculate  of feedback effects by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)
    # Extract normalized scores (divided by the column's mean) to get a feedback index for each construct polarity
    feedback <- apply (feedback, MARGIN = 1, FUN = mean)
  }

  # Define the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (feedback , consequence)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["constructs"]]
  # Show row names (titles)
  rownames(result) <- poles

  return(result)
}

# cq.fb2 ---------------------------------------------------------------

#' Calculation of consequence and feedback effects for each construct of the
#' system
#'
#' @description This function calculates the numbers of consequences and feedbacks
#' @description This function calculates the intensity of the effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp
#' @param std
#' @param std A logical value indicating wether to standarize the scores (TRUE) or (FALSE)
#'
#' @return Returns a dataframe with numbers of consequences and feedbacks
#' @return Returns a dataframe with the intensity of the effects of consequences and feedbacks
#'
#'
#'
#' @export
#'
#' @examples
#' cq.fb2_index (example.wimp)
#'

cq.fb2_index <- function(wimp) {

  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights

  # Extracting scores of positive consequences as measured when comparing to the hypothetical self
  cons_pos <- numeric(nrow(wmatrix))

  for (i in 1:nrow(wmatrix)) {
    fila <- wmatrix[i, ]
    cons_pos[i] <- sum(fila[fila > 0], na.rm = TRUE)
  }

  # Extracting scores of negative consequences as measured when comparing to the hypothetical self
  cons_neg <- numeric(nrow(wmatrix))

  for (i in 1:nrow(wmatrix)) {
    fila <- wmatrix[i, ]
    cons_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)
  }

    # Calculate scores of feedback effects by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)

    # Extract scores of positive feedback index for each construct polarity
    fb_pos <- numeric(nrow(feedback))

    for (i in 1:nrow(feedback)) {
      fila <- feedback[i, ]
      fb_pos[i] <- sum(fila[fila > 0], na.rm = TRUE)
    }

    # Extract scores of negative feedback index for each construct polarity
    fb_neg <- numeric(nrow(feedback))

    for (i in 1:nrow(feedback)) {
      fila <- feedback[i, ]
      fb_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)
    }


  # Define the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (cons_neg, fb_neg, cons_pos, fb_pos)
  # Setting names to the poles of each construct
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
  # Calling function cq.fb_index
  result <- cq.fb_index (wimp, std = std)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["constructs"]]
  # Show row names (titles)
  rownames(result) <- poles
  # Create a scatter plot to visualize the effects of cq and fb
  fig <- plot_ly(data = result, x = ~consequence, y = ~feedback, color = ~poles)
  return(fig)
}
