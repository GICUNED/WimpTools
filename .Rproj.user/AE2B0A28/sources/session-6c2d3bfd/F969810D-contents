## CHANGE IMPLICATIONS FUNCTIONS ##

# IF Index function---------------------------------------------------------------

#' Impact And Feedback Index -- if_index()
#'
#' @description This function calculates the intensity of the of impact in the
#' system and feedback from the system for each construct.
#'
#' @param wimp An object of class 'wimp', which contains a Weigthed Implication
#'             Grid.
#' @param std A character string indicating the type of standardization to apply.
#'        Options are 'none' (no standardization), 'vertex' (standardize by the maximum total degree
#'        of the constructs), 'edges' (standardize by the total number of edges) and 'adjacent'
#'        (standarized by the number of the adjacent vertices). Default is 'adjacent'.
#'
#' @author Maite Benitez Santos, Guillermo Calleja Garate, Alejandro Sanfeliciano
#'
#' @return Returns a dataframe showing the intensity of the impact in the
#' system and feedback from the system for each construct (negative, postive and
#'  global).
#'
#' @export
#'
#' @examples
#'
#' if_index (example.wimp)
#'

if_index <- function(wimp, std = "adjacent") {

  # Aligning the wimp object to ideal
  wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

  # Extract weight matrix from the wimp
  wmatrix <- wimp$scores$weights
  amatrix <- wmatrix
  amatrix[amatrix != 0] <- 1
  ncol <- ncol(wmatrix)

  # Define standarization coefficient

  if(std == "vertex"){std.coef <- dim(wmatrix)[1] - 1}
  if(std == "edges"){std.coef <- length(wmatrix[wmatrix != 0])}
  if(std == "adjacent"){std.coef <- rowSums(amatrix)}

  self.vector <- wimp$self$standarized
  ideal.vector <- wimp$ideal$standarized

  self.matrix <- matrix(rep(self.vector, ncol), nrow = length(self.vector), ncol = ncol, byrow = TRUE)
  ideal.matrix <- matrix(rep(ideal.vector, ncol), nrow = length(ideal.vector), ncol = ncol, byrow = TRUE)

  # Impact calculation
  imatrix <- wmatrix
  imatrix[,.which.dilemmatics(wimp)] <- 0

  #imatrix <- mapply(.impact, w = wmatrix, s=self.matrix , i=ideal.matrix)
  #imatrix <- matrix(imatrix , ncol = ncol)

  imatrix.positive <- imatrix
  imatrix.positive[imatrix.positive < 0] <- 0

  imatrix.negative <- imatrix
  imatrix.negative[imatrix.positive > 0] <- 0

  negative.impact <- apply(imatrix.negative, MARGIN = 1, FUN = sum)
  positive.impact <- apply(imatrix.positive, MARGIN = 1, FUN = sum)
  global.impact <- apply(imatrix, MARGIN = 1, FUN = sum)

  if(std != "none"){
    negative.impact <- negative.impact / std.coef
    positive.impact <- positive.impact / std.coef
    global.impact <- global.impact / std.coef
  }

  # Feedback calculation

  fmatrix <- wmatrix * t(wmatrix)

  fmatrix.positive <- fmatrix
  fmatrix.positive[fmatrix.positive < 0] <- 0

  fmatrix.negative <- fmatrix
  fmatrix.negative[fmatrix.positive > 0] <- 0

  negative.feedback <- apply(fmatrix.negative, MARGIN = 1, FUN = sum)
  positive.feedback  <- apply(fmatrix.positive, MARGIN = 1, FUN = sum)
  global.feedback  <- apply(fmatrix, MARGIN = 1, FUN = sum)

  if(std != "none"){
    negative.feedback <- negative.feedback / std.coef
    positive.feedback <- positive.feedback / std.coef
    global.feedback <- global.feedback / std.coef
  }

  # Defining data.frame
  result <- data.frame(abs(negative.impact),positive.impact, global.impact,
                        abs(negative.feedback),positive.feedback,
                        global.feedback)

  # Set up col and row names
    colnames(result) <- c("Negative Impact", "Positive Impact", "Global Impact",
                          "Negative Feedback", "Positive Feedback",
                          "Global Feedback")
    rownames(result) <- wimp$constructs$constructs

  # Return function
    return(result)
}

# IF Plot Function ----------------------------------------------------------------

#' IF Index Plot  -- if_plot()
#'
#' @description This function creates a scatter plot to show the results of the
#'              \code{\link{if_index}} function.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param show A character string indicating the constructs to be displayed. 'all' will display
#'        all constructus, 'dil' will only display dilemmatic constructs and 'nodil' will
#'        exclude dilemmatic constructs. Default is 'all'.
#' @param text.size Scalar that modifies the text size. Default is 1.
#' @param center Establishes the centre of the frame. Use "data" to set the data
#'        to be framed and "origin" to set the origin to be in the centre. the default
#'        is "data".
#' @param ... additional arguments are passed from \code{\link{if_index}}
#'        function.
#'
#' @author Maite Benitez Santos, Guillermo Calleja Garate and Alejandro Sanfeliciano
#'
#' @return returns a interactive scatter plot made with Plotly.
#'
#' @export
#'
#' @import plotly
#'
#' @examples
#'
#' if_plot (example.wimp)
#'

if_plot <- function(wimp, show = "all", center = "data", text.size = 1, ...) {

    # Align wimpgrid towards ideal
    wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

    # Extract important info
    self.poles <- wimp$constructs$self.poles
    right.poles <- wimp$constructs$right.poles
    constructs <- wimp$constructs$constructs
    dil <- .which.dilemmatics(wimp)
    construct.color <- .construct.colors(wimp, "red/green")

    # Set up data.frame
    df <- if_index(wimp, ...)[c(3,6)]
    df <- data.frame (df, right.poles, constructs, self.poles, construct.color)


    # Row and col names for data.frame
    names(df) <- c("I", "FB", "poles", "construct", "self", "color")
    rownames(df) <- right.poles

    # Show options
    if(show == "nodil"){
      df <- df[-dil,]
    }
    if(show == "dil"){
      df <- df[dil,]
    }

    # Plotting
    # Axis set up
    if(center == "data"){
    irange.min <- min(df[1]) - 0.15 * max(abs(df[1]))
    irange.max <- max(df[1]) + 0.15 * max(abs(df[1]))

    frange.min <- min(df[2]) - 0.15 * max(abs(df[2]))
    frange.max <- max(df[2]) + 0.15 * max(abs(df[2]))
    }

    if(center == "origin"){

      irange.min <- -(max(abs(df[1])) + 0.15 * max(abs(df[1])))
      irange.max <- max(abs(df[1])) + 0.15 * max(abs(df[1]))

      frange.min <- -(max(abs(df[2])) + 0.15 * max(abs(df[2])))
      frange.max <- max(abs(df[2])) + 0.15 * max(abs(df[2]))
    }

    # Shapes
    shapes <- list(
      list(type = "rect",
           fillcolor = "palegreen", line = list(color = "palegreen"), opacity = 0.3, layer="below",
           x0 = 0, x1 = 10000,
           y0 = 0, y1 = 10000),

      list(type = "rect",
           fillcolor = "#eb636b", line = list(color = "#eb636b"), opacity = 0.3, layer="below",
           x0 = 0, x1 = -10000,
           y0 = 0, y1 = -10000)
      ,
      list(type = "rect",
           fillcolor = "#ffe65d", line = list(color = "#ffe65d"), opacity = 0.3, layer="below",
           x0 = 0, x1 = -10000,
           y0 = 0, y1 = 10000),

      list(type = "rect",
           fillcolor = "#ffe65d", line = list(color = "#ffe65d"), opacity = 0.3, layer="below",
           x0 = 0, x1 = 10000,
           y0 = 0, y1 = -10000)
    )

    # Scatter plot
    fig <- plot_ly(
      data = df,
      x = ~I,
      y = ~FB
    ) %>% add_annotations(
      data = df,
      x = ~I,
      y = ~FB,
      text = ~poles,
      hoverinfo = 'text',
      font = list(size = 15 * text.size),
      showarrow = FALSE,
      xanchor = 'center',
      yanchor = 'bottom',
      yshift = 5
    ) %>% add_markers(
      data = df,
      x = ~I,
      y = ~FB,
      marker = list(color = ~color, size = 7, line = list(color = 'black', width = 1)),
      text = ~paste('<B>',construct,'</B>' , '\nSelf:', self, '\nI:', round(I,digits = 2), '\nF:', round(FB,digits = 2)),
      hoverinfo = 'text'
    ) %>% layout(
             xaxis = list(title = "IMPACT",
                          range = c(irange.min,irange.max),
                          gridcolor ="white",
                          gridwidth= 0.5,
                          zeroline = TRUE,
                          zerolinecolor = "black",
                          zerolinewidth = 2),
             yaxis = list(title = "FEEDBACK",
                          range = c(frange.min,frange.max),
                          gridcolor ="white",
                          gridwidth= 0.5,
                          zeroline = TRUE,
                          zerolinecolor = "black",
                          zerolinewidth = 2),
             showlegend = FALSE,
             shapes = shapes
             )

    # Return function
    return(fig)
}

# Impact and Feedback Bar Chart Function ----------------------------------------------------------------


#' Impact and Feedback Bar Chart  -- if_barchart()
#'
#' @description This function creates a bar chart to show the results of the
#'              \code{\link{if_index}} function.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param show A character string indicating the constructs to be displayed. 'all' will display
#'        all constructus, 'dil' will only display dilemmatic constructs and 'nodil' will
#'        exclude dilemmatic constructs. Default is 'all'.
#' @param ... additional arguments are passed from \code{\link{if_index}}
#'        function.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return returns a interactive bar chart made with plotly.
#'
#' @export
#'
#' @import plotly
#'
#' @examples
#'
#' if_barchart(example.wimp)
#'

if_barchart <- function(wimp, show = "all",...){

  # Align wimp
  wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

  #Extract important info

  self.poles <- wimp$constructs$self.poles
  rigth.poles <- wimp$constructs$right.poles
  constructs <- wimp$constructs$constructs

  dim <- length(constructs)
  dil <- .which.dilemmatics(wimp)

  color.line.p <- rep("#5ce75c",dim)
  if(show == "all"){color.line.p[dil] <- "#FFD97D"}

  color.line.n <- rep("#d13b43",dim)
  if(show == "all"){color.line.n[dil] <- "#FFD97D"}

  width.line <- rep(1,dim)
  if(show == "all"){width.line[dil] <- 1}

  pattern <- rep(0,dim)
  pattern[dil] <- 1
  # Set up data frame for further plotting
  df <- if_index(wimp,...)[c(1,2,4,5)]
  df <- data.frame(df,df[,1] + df[,2],rigth.poles,self.poles,constructs,color.line.n,color.line.p,width.line)
  df[,1] <- -df[,1]
  df[,3] <- -df[,3]

  # Row and col names of data frame
  names(df) <- c("NI","PI","NF","PF","global","right.poles","self.poles","construct","colorn","colorp","widthline")
  rownames(df) <- rigth.poles

  # Show options
  if(show == "nodil"){
    df <- df[-dil,]
  }
  if(show == "dil"){
    df <- df[dil,]
  }

  # Plotting
  # Important info to plot
  range <- max(abs(df[c(1,2,3,4)])) + 0.15 * max(abs(df[c(1,2,3,4)]))

  # Plot fig1
  fig1 <- plot_ly(
    data = df,
    x = ~PI,
    y = ~reorder(right.poles, global),
    type = "bar",
    orientation = "h",
    marker = list(color = "#AAF683", line = list(color = ~colorp, width = ~widthline), pattern = list(shape = ~ifelse(pattern == 1, "/", ""),fillmode="overlay", fgcolor="#FFEE7D", size=20)),
    hovertext = ~paste('<B>', construct,'</B>', '\nSelf:', self.poles, '\nPositive Impact:', round(PI, 2)),
    hoverinfo = "text"
  ) %>% add_trace(
    x = ~NI,
    name = 'Impact',
    marker = list(color = "#EE6055", line = list(color = ~colorn, width = ~widthline), pattern = list(shape = ~ifelse(pattern == 1, "/", ""),fillmode="overlay", fgcolor="#FFEE7D", size=20)),
    hovertext = ~paste('<B>', construct,'</B>','\nSelf:', self.poles, '\nNegative Impact:', round(NI, 2)),
    hoverinfo = "text"
  ) %>% layout(
    barmode = "overlay",
    bargap = 0.08,
    xaxis = list(title = "IMPACT", range = c(-range, range), showline = TRUE),
    yaxis = list(title = "", showgrid = TRUE, showline = TRUE),
    showlegend = FALSE
  )

  # Plot fig2
  fig2 <- plot_ly(
      data = df,
      x = ~PF,
      y = ~reorder(right.poles, global),
      type = "bar",
      orientation = "h",
      name = "Feedback",
      marker = list(color = "#AAF683", line = list(color = ~colorp, width = ~widthline), pattern = list(shape = ~ifelse(pattern == 1, "/", ""),fillmode="overlay", fgcolor="#FFEE7D", size=20)),
      hovertext = ~paste('<B>', construct, '</B>', '\nSelf:', self.poles, '\nPositive Feedback:', round(PF, 2)),
      hoverinfo = "text"
  ) %>% add_trace(
      x = ~NF,
      name = 'Feedback',
      marker = list(color = "#EE6055", line = list(color = ~colorn, width = ~widthline), pattern = list(shape = ~ifelse(pattern == 1, "/", ""),fillmode="overlay", fgcolor="#FFEE7D", size=20)),
      hovertext = ~paste('<B>', construct,'</B>', '\nSelf:', self.poles, '\nNegative Feedback:', round(NF, 2)),
      hoverinfo = "text"
  ) %>% layout(
      barmode = "overlay",
      bargap = 0.08,
      xaxis = list(title = "FEEDBACK", range = c(-range, range), showline = TRUE),
      yaxis = list(title = "", showgrid = TRUE, showline = TRUE, showticklabels = TRUE, side = "right"),
      showlegend = FALSE
  )

  # Combine fig1 and f2 in a subplot
  fig <- subplot(
      fig1, fig2, margin = 0.005
  ) %>% layout(
      xaxis = list(title = "IMPACT",
                   gridcolor = "white",
                   gridwidth = 3),
      yaxis = list(title = "",
                   gridcolor = "white",
                   gridwidth = 3),
      xaxis2 = list(title = "FEEDBACK",
                   gridcolor = "white",
                   gridwidth = 3),
      yaxis2 = list(title = "", showticklabels = TRUE, side = "right", overlaying = "y",
                   gridcolor = "white",
                   gridwidth = 3),
      plot_bgcolor = "#FFFAED"
  )

  # Return function
  return(fig)
  }

