## RESISTANCE TO CHANGE FUNCTIONS ##



# cq.fb ---------------------------------------------------------------



#' Calculation of consequence and feedback effects for each construct of the system in three different ways:
#' Mode 1: Calculation the intensity of the effects of consequences and feedbacks
#' Mode 2: Calculation of consequence and feedback negative effects for each construct of the system
#' Mode 3: Calculation of consequence and feedback of both positive and negative effects for each construct of the system
#'
#' @description This function calculates the intensity of the effects of consequences and feedback
#' for each construct to analyze its resistance to change in three different ways
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and its associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#' @param mode A value with three options: 1 = Sum of consequences and feedback; 2 = Sum of negative consequences and feedback;
#' 3 = Sum of positive and negative consequences and feedback
#' @return Returns a dataframe showing the intensity of the effects caused by consequences and feedback interactions according to the previously selected mode
#'
#' @export
#'
#' @examples
#' cq.fb_index (example.wimp)
#'



cq.fb_index <- function(wimp, mode = 1, std = FALSE) {



  # Aligning the wimp object so that ideal situations get close to the right pole

  wimp <- .align.wimp(wimp)

  # Creating a weight matrix from the wimp to be analysed

  wmatrix <- wimp $scores $weights

  # Setting up modes and standarization

  if (mode == 1) {

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

  else if (mode == 2) {

  # Calculation of consequence and feedback negative effects for each construct of the system

    if (std == FALSE) {

  # Extracting scores of negative consequences as measured when comparing to the hypothetical self

      cons_neg <- numeric(nrow(wmatrix))

      for (i in 1:nrow(wmatrix)) {

        fila <- wmatrix[i, ]

        cons_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)

        if (is.nan(cons_neg[i])){cons_neg[i] <-0}
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

        if (is.nan(cons_neg[i])){cons_neg[i] <-0}
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

  } else {

  # Calculation of consequence and feedback of positive and negative effects for each construct of the system

    if (std == FALSE) {

  # Extracting scores of positive consequences as measured when comparing to the hypothetical self

      cons_pos <- numeric(nrow(wmatrix))

      for (i in 1:nrow(wmatrix)) {

        fila <- wmatrix[i, ]

        cons_pos[i] <- sum(fila[fila > 0], na.rm = TRUE)

        if (is.nan(cons_pos[i])){cons_pos[i] <-0}
      }

  # Extracting scores of negative consequences as measured when comparing to the hypothetical self

      cons_neg <- numeric(nrow(wmatrix))

      for (i in 1:nrow(wmatrix)) {

        fila <- wmatrix[i, ]

        cons_neg[i] <- sum(fila[fila < 0], na.rm = TRUE)

        if (is.nan(cons_neg[i])){cons_neg[i] <-0}
      }

  # Calculating feedback scores by multiplying the weight matrix by its own transposed version

      feedback <- wmatrix * t (wmatrix)

  # Extracting scores of positive feedback index for each construct polarity

      fb_pos <- numeric(nrow(feedback))

      for (i in 1:nrow(feedback)) {

        fila <- feedback[i, ]

        fb_pos[i] <- sum(fila[fila > 0], na.rm = TRUE)

        if (is.nan(fb_pos[i])){fb_pos[i] <-0}
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

        if (is.nan(cons_pos[i])){cons_pos[i] <-0}
      }

  # Extracting scores of negative consequences as measured when comparing to the hypothetical self

      cons_neg <- numeric(nrow(wmatrix))

      for (i in 1:nrow(wmatrix)) {

        fila <- wmatrix[i, ]

        cons_neg[i] <- mean(fila[fila < 0], na.rm = TRUE)

        if (is.nan(cons_neg[i])){cons_neg[i] <-0}
      }

  # Calculating feedback scores by multiplying the weight matrix by its own transposed version

      feedback <- wmatrix * t (wmatrix)

  # Extracting positive feedback score indices for each construct polarity

      fb_pos <- numeric(nrow(feedback))

      for (i in 1:nrow(feedback)) {

        fila <- feedback[i, ]

        fb_pos[i] <- mean(fila[fila > 0], na.rm = TRUE)

        if (is.nan(fb_pos[i])){fb_pos[i] <-0}
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
}






# rtc_plot ----------------------------------------------------------------
#' Visualization of the effects of Consequences and Feedback using two different ways
#'
#' @description This function creates a scatter plot to show the effects of the consequences and feedbacks of each construct to analyze their resistance to change
#' keeping in mind that mode 2 only shows the negative effects
#'
#' @param wimp An object of class 'wimp', which contains an implication grid and associated constructs.
#' @param std A logical value indicating whether to standarize the scores (TRUE) or (FALSE)
#' @param mode Value with two options: 1 = shows the sum of effects of consequences and feedback on each construct of the system;
#' 2= shows negative effects only
#' @return returns a scatter plot showing the effects of consequences and feedback on each construct of the system
#' @export
#' @import plotly
#' @examples
#' rtc_plot (example.wimp)

rtc_plot <- function(wimp, mode=1, ...) {

  if (mode == 1) {
    # Calling function cq.fb_index

    result <- cq.fb_index (wimp, mode=mode, ...)
    construct_color <- construct_colors(wimp, "red/green")

    # Defining the data frame needed for numeric visualization of effects and further plotting

    result <- data.frame (result, construct_color)

    # Setting names to the poles of each construct

    poles <- wimp[["constructs"]][["self.poles"]]
    bothpoles <- wimp[["constructs"]][["constructs"]]

    # Showing row names (titles)

    rownames(result) <- poles
    rango <- max(result[1:2]) + 0.05

    # Creating a scatter plot to visualize the effects of consequences and feedback

    fig <- plot_ly(data = result [1:2], x = ~consequence, y = ~feedback, color= ~construct_color)

    # Adding labels to graph points and axes

    fig <- fig %>%

      add_trace(type = "scatter", mode = "markers", text = bothpoles, showlegend = FALSE) %>%

      add_annotations(data = result,text = poles,

                      hoverinfo = 'text',

                      font = list(size = 12),

                      showarrow = TRUE, xanchor = 'center', yanchor = 'bottom',

                      yshift = 5)  %>%

      layout(title = "CONSEQUENCES AND FEEDBACK EFFECTS",

             xaxis = list(title = "Consequences", range=c(-rango,rango)),

             yaxis = list(title = "Feedback", range=c(-rango,rango)) )

    return(fig)
  }

  if (mode == 2) {

  # Visualization of the negative effects of Consequences and Feedback
  # Calling function cq.fb_index

    result <- cq.fb_index (wimp, mode=mode, ...)
    construct_color <- construct_colors(wimp, "red/green")
    result <- data.frame (result, construct_color)

  # Setting names to the poles of each construct

    poles <- wimp[["constructs"]][["self.poles"]]
    bothpoles <- wimp[["constructs"]][["constructs"]]

  # Showing row names (titles)

    rownames(result) <- poles

  # Creating a scatter plot to visualize negative effects of consequences and feedback

    fig <- plot_ly(data = result[1:2], x = ~cons_neg, y = ~fb_neg, color = construct_color)

  # Adding labels to graph points and axes

    fig <- fig %>%

      add_trace(type = "scatter", mode = "markers", text = bothpoles, showlegend = FALSE) %>%
      add_annotations(data = result[1:2],text = poles,

                      hoverinfo = 'text',

                      font = list(size = 12, color = construct_color),

                      showarrow = TRUE, xanchor = 'center', yanchor = 'bottom',

                      yshift = 5)  %>%

      layout(title = "NEGATIVE EFFECTS OF CONSEQUENCES AND FEEDBACK",

             xaxis = list(title = "Negative Consequences" ),

             yaxis = list(title = "Negative Feedback"))

  # Returning the plot

    return(fig)
  }
}

