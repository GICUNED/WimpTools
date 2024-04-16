## RESISTANCE TO CHANGE FUNCTIONS ##

# cq.fb ---------------------------------------------------------------

#' Calculation of consequence and feedback effects for each construct of the
#' system
#'
#' @description This function calculates the intensity of the effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#'
#' @return Returns a dataframe with the intensity of the effects caused by consequences and feedback interactions
#'
#' @export
#'
#' @examples
#' cq.fb_index (example.wimp)
#'

cq.fb_index <- function(wimp, std = FALSE) {

  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights

  if (std == FALSE) {
    # Extracting scores of consequences as measured when comparing to the hypothetical self
    consequence <- apply (wmatrix, MARGIN = 1, FUN = sum)
    # Calculating feedback scores by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)
    # Extracting feedback scores to get a feedback index for each construct polarity
    feedback <- apply (feedback, MARGIN = 1, FUN = sum)

  } else {
    # Extracting normalized scores (divided by the column's mean) of consequences as measured when comparing to the hypothetical self
    consequence <- apply (wmatrix, MARGIN = 1, FUN = mean)
    # Calculating  of feedback effects by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)
    # Extracting normalized scores (divided by the column's mean) to get a feedback index for each construct polarity
    feedback <- apply (feedback, MARGIN = 1, FUN = mean)
  }

  # Defining the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (feedback , consequence)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["self.poles"]]
  # Showing row names (titles)
  rownames(result) <- poles

  return(result)
}

# cq.fb2 ---------------------------------------------------------------

#' Calculation of consequence and feedback of positive and negative effects for each construct of the
#' system
#'
#'
#' @description This function calculates the intensity of positive and negative effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#'
#'
#' @return Returns a dataframe with the intensity of positive and negative effects caused by consequences and feedback interactions
#'
#'
#'
#' @export
#'
#' @examples
#' cq.fb2_index (example.wimp)
#'

cq.fb2_index <- function(wimp, std = FALSE) {

  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights

  if (std == FALSE) {

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

  # Calculating feedback scores by multiplying the weight matrix by its own transposed version
  feedback <- wmatrix * t (wmatrix)

  # Extracting scores of positive feedback index for each construct polarity
  fb_pos <- numeric(nrow(feedback))

  for (i in 1:nrow(feedback)) {
    fila <- feedback[i, ]
    fb_pos[i] <- sum(fila[fila > 0], na.rm = TRUE)
  }

  # Extracting scores of negative feedback indices for each construct polarity
  fb_neg <- numeric(nrow(feedback))

  for (i in 1:nrow(feedback)) {
    fila <- feedback[i, ]
    fb_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)
  }

  } else {

    # Extracting scores of positive consequences as measured when comparing to the hypothetical self
    cons_pos <- numeric(nrow(wmatrix))

    for (i in 1:nrow(wmatrix)) {
      fila <- wmatrix[i, ]
      cons_pos[i] <- mean(fila[fila > 0], na.rm = TRUE)
    }

    # Extracting scores of negative consequences as measured when comparing to the hypothetical self
    cons_neg <- numeric(nrow(wmatrix))

    for (i in 1:nrow(wmatrix)) {
      fila <- wmatrix[i, ]
      cons_neg[i] <- mean(fila[fila < 0], na.rm = TRUE)
    }

    # Calculating feedback scores by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)

    # Extracting positive feedback score indices for each construct polarity
    fb_pos <- numeric(nrow(feedback))

    for (i in 1:nrow(feedback)) {
      fila <- feedback[i, ]
      fb_pos[i] <- mean(fila[fila > 0], na.rm = TRUE)
    }

    # Extracting negative feedback score indices for each construct polarity
    fb_neg <- numeric(nrow(feedback))

    for (i in 1:nrow(feedback)) {
      fila <- feedback[i, ]
      fb_neg[i] <- mean(fila[fila < 0], na.rm = TRUE)
    }
  }

  # Defining the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (cons_neg, fb_neg, cons_pos, fb_pos)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["self.poles"]]
  # Showing row names (titles)
  rownames(result) <- poles
  return(result)
}

# cq.fb3 ---------------------------------------------------------------

#' Calculation of consequence and feedback negative effects for each construct of the
#' system
#'
#' @description This function calculates the intensity of negative effects of consequences and feedbacks
#' for each construct to analyze its resistance to change
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#'
#' @return Returns a dataframe with the intensity of negative effects caused by consequences and feedback interactions
#'
#'
#'
#' @export
#'
#' @examples
#' cq.fb3_index (example.wimp)
#'


cq.fb3_index <- function(wimp, std = FALSE) {

  # Aligning the wimp object so that ideal situations get close to the right pole
  wimp <- .align.wimp(wimp)
  # Creating a weight matrix from the wimp to be analysed
  wmatrix <- wimp $scores $weights

  if (std == FALSE) {
  # Extracting scores of negative consequences as measured when comparing to the hypothetical self
  cons_neg <- numeric(nrow(wmatrix))

  for (i in 1:nrow(wmatrix)) {
    fila <- wmatrix[i, ]
    cons_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)
  }

  # Calculating feedback scores by multiplying the weight matrix by its own transposed version
  feedback <- wmatrix * t (wmatrix)


  # Extracting negative feedback score indices for each construct polarity
  fb_neg <- numeric(nrow(feedback))

  for (i in 1:nrow(feedback)) {
    fila <- feedback[i, ]
    fb_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)
  }

  } else {

    cons_neg <- numeric(nrow(wmatrix))
    for (i in 1:nrow(wmatrix)) {
      fila <- wmatrix[i, ]
      cons_neg[i] <- mean(fila[fila < 0], na.rm = TRUE)
    }

    # Calculating feedback scores by multiplying the weight matrix by its own transposed version
    feedback <- wmatrix * t (wmatrix)


    # Extracting negative feedback score indices for each construct polarity
    fb_neg <- numeric(nrow(feedback))

    for (i in 1:nrow(feedback)) {
      fila <- feedback[i, ]
      fb_neg[i] <- mean(fila[fila < 0], na.rm = TRUE)
    }
  }

  # Defining the data frame needed for numeric visualization of effects and further plotting
  result <- data.frame (cons_neg, fb_neg)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["self.poles"]]
  # Showing row names (titles)
  rownames(result) <- poles
  return(result)
}



# rtc_plot ----------------------------------------------------------------


#' Visualization of the effects of Consequences and Feedback
#'
#' @description This function creates a scatter plot to show the effects of the consequences and feedbacks of each construct to analyze their resistance to change
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
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
  poles <- wimp[["constructs"]][["self.poles"]]
  # Showing row names (titles)
  rownames(result) <- poles
  # Creating a scatter plot to visualize the effects of consequences and feedback
  fig <- plot_ly(data = result, x = ~consequence, y = ~feedback, color = ~poles)

  # Adding labels to graph points and axes
  fig <- fig %>%
    add_trace(type = "scatter", mode = "markers", text = poles, showlegend = FALSE) %>%
    layout(title = "Consequences and Feedback Effects",
           xaxis = list(title = "Consequences"),
           yaxis = list(title = "Feedback"))

  return(fig)
}




# rtc2_plot ----------------------------------------------------------------


#' Visualization of the negative effects of Consequences and Feedback
#'
#' @description This function creates a scatter plot to show the negative effects of the consequences and feedbacks of each construct to analyze their resistance to change
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#'
#' @return returns a scatter plot showing the negative effects of consequences and feedback on each construct of the system
#' @export
#' @import plotly
#' @examples
#' rtc2_plot (example.wimp)

rtc2_plot <- function(wimp, std = FALSE) {
  # Calling function cq.fb2_index
  result <- cq.fb3_index (wimp, std = std)
  # Setting names to the poles of each construct
  poles <- wimp[["constructs"]][["self.poles"]]
  # Showing row names (titles)
  rownames(result) <- poles

  # Creating a scatter plot to visualize negative effects of consequences and feedback
  fig <- plot_ly(data = result, x = ~cons_neg, y = ~fb_neg, color = ~poles)

  # Adding labels to graph points and axes
  fig <- fig %>%
    add_trace(type = "scatter", mode = "markers", text = poles, showlegend = FALSE) %>%
    layout(title = "Negative Effects of Consequences and Feedback",
           xaxis = list(title = "Negative Consequences"),
           yaxis = list(title = "Negative Feedback"))

  # Returning the plot
  return(fig)


}
