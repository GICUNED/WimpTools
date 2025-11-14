## ADJUSTMENT FUNCTIONS ##

# Construct proportions index  ---------------------------------------------------

#' Frencuency and proportions of constructs -- construct_index()
#'
#' @description This function calculates frequency and proportion of
#'              congruents, discrepants, dilemmatics and undefined constructs.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A matrix with the frequency and proportion of congruents, discrepants
#'        , dilemmatics and undefined constructs.
#'
#' @export
#'
#' @examples
#'
#' construct_index(example.wimp)
#'

construct_index <- function(wimp){

indices <- .wimp_get_construct_indices(wimp)
n.congruents <- length(indices$congruents)
n.discrepants <- length(indices$discrepants)
n.dilemmatics <- length(indices$dilemmatics)
n.undefined <- length(indices$undefined)

n <- .wimp_n_constructs(wimp)

congruents <- c(n.congruents, n.congruents / n)
discrepants <- c(n.discrepants, n.discrepants / n)
dilemmatics <- c(n.dilemmatics, n.dilemmatics / n)
undefined <- c(n.undefined, n.undefined / n)

result <- rbind(congruents,discrepants,dilemmatics,undefined)

rownames(result) <- c("Congruents","Discrepants","Dilemmatics","Undefined")
colnames(result) <- c("Frequency","Proportion")

return(result)
}


# Self Correlations ---------------------------------------------------

#' Correlations between Self and Hypothetical scenarios -- self_index()
#'
#' @description this function Calculates the global and specific adjustment
#'              indices of the self for each hypothetical scenario in the WimpGrid.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param method The correlation method to use. All methods of the \code{\link{cor}} function
#'        are allowed and "ssi" for SSI Index. Default is "ssi".
#' @param rc Use Cohen's rc which is invariant to construct reflection. Default is TRUE.
#' @param alpha Alpha value for SSI calculation. Default is .5.
#' @param beta Beta value for SSI calculation. Default is .5.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A list with global adjustment self indices and specific indices for each construct.
#'
#' @export
#'
#' @examples
#'
#' self_index(example.wimp)
#'

self_index <- function(wimp, method = "ssi", rc = TRUE, alpha = .5, beta = .5){

  result <- list()

  indices <- .wimp_get_construct_indices(wimp)
  congruence <- rep("Undefined", .wimp_n_constructs(wimp))
  congruence[indices$congruents] <- "Discrepant"
  congruence[indices$undefined] <- "Congruent"
  congruence[indices$discrepants] <- "Congruent"
  congruence[indices$dilemmatics] <- "Undefined"

  hypo.matrix <- .hypo.matrix(wimp)
  ncol <- ncol(hypo.matrix)
  hypo.names <- colnames(hypo.matrix)[-c(1,ncol)]
  hypo.names <- paste("Totally", hypo.names, sep = " ")

  rc.text <- "no rc"
  if(rc){
    hypo.matrix <- rbind(hypo.matrix, -hypo.matrix)
    rc.text <- "rc"
    }

  if(method == "ssi"){

    self.vector <- hypo.matrix[,1]
    ideal.vector <- hypo.matrix[,ncol]

    discrepants <- indices$discrepants
    congruents <- indices$congruents

    self.cor <- apply(hypo.matrix[,-c(1, ncol(hypo.matrix))], 2, function(col) .sim_index(self.vector, col, alpha = alpha, beta = beta))
    ideal.cor <- apply(hypo.matrix[,-c(1, ncol(hypo.matrix))], 2, function(col) .sim_index(ideal.vector, col, alpha = alpha, beta = beta))

    self.ideal.cor <- .sim_index(self.vector,ideal.vector, alpha = alpha, beta = beta)
    self.hypo.cor <- mean(self.cor)
    ideal.hypo.cor <- mean(ideal.cor)

    ideal.hypo.congruents.cor <- mean(ideal.cor[congruents])
    ideal.hypo.discrepants.cor <- mean(ideal.cor[discrepants])

  }
  if(!( method == "ssi" | method == "cos" )){
    self.vector <- hypo.matrix[,1]
    ideal.vector <- hypo.matrix[,ncol]

    discrepants <- wimp$constructs$discrepants
    congruents <- wimp$constructs$congruents

    self.cor <- cor(self.vector,hypo.matrix[,-c(1,ncol)], method = method)
    ideal.cor <- cor(ideal.vector,hypo.matrix[,-c(1,ncol)], method = method)

    self.ideal.cor <- cor(self.vector,ideal.vector, method = method)
    self.hypo.cor <- mean(self.cor)
    ideal.hypo.cor <- mean(ideal.cor)

    ideal.hypo.congruents.cor <- mean(ideal.cor[congruents])
    ideal.hypo.discrepants.cor <- mean(ideal.cor[discrepants])
  }

  global <- mean(cor(ideal.vector,hypo.matrix[,-ncol]))

  df.global <- data.frame(self.ideal.cor,self.hypo.cor,ideal.hypo.cor,ideal.hypo.congruents.cor,ideal.hypo.discrepants.cor)
  names(df.global) <- c("Self/Ideal", "Self/Hypo", "Ideal/Hypo", "Ideal/Discrepant", "Ideal/Congruent")

  df.construct <- data.frame(
    Hypothetical_Scenario = hypo.names,
    Congruence_Scenario = congruence,
    Self_Similarity = round(as.numeric(self.cor), 4),
    Ideal_Similarity = round(as.numeric(ideal.cor), 4),
    stringsAsFactors = FALSE
  )

  names(df.construct) <- c("Hypothetical Scenario","Congruence Scenario","Self Similarity", "Ideal Similarity")

  result$global <- df.global
  result$construct <- df.construct
  result$method <- c(method, rc.text)

  return(result)
}

# Adjustment Radar Chart ---------------------------------------------------

#' Adjustment Radar Chart -- adj_plot()
#'
#' @description This function creates a radar chart showing the value of the
#'              self for each construct and its adjustment with respect to the ideal.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by  the \code{\link{importwimp}} function.
#'
#' @return A Plotly radar polar plot.
#'
#' @author Alejandro Sanfeliciano
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#' adj_plot(example.wimp)
#'

adj_plot <- function(wimp){

  wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

  self <- .wimp_get_self(wimp)
  self <- c(self,self[1])

  ideal <- .wimp_get_ideal(wimp)
  ideal <- c(ideal,ideal[1])

  r.poles <- .wimp_get_right_poles(wimp)
  l.poles <- .wimp_get_left_poles(wimp)
  poles <- paste(r.poles," (",l.poles,")", sep="")
  poles <- c(poles,poles[1])

  construct<- .wimp_get_construct_names(wimp)
  construct <- c(construct,construct[1])

  colors <- .construct.colors(wimp, mode = "red/green")[,1]
  colors <- c(colors,colors[1])

  plot <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  )
  plot <- plot %>%
    add_trace(
      mode = "lines",
      r = 0,
      theta = poles,
      fill = "none",
      line = list(color = "#444444", width = 1.5, shape = 'spline', smoothing = 1),
      name = 'Pole Threshold',
      hoverinfo = 'none'
    )
  plot <- plot %>%
    add_trace(
      mode = "lines",
      r = ideal,
      theta = poles,
      fill = "none",
      line = list(color = "darkgreen", width = 3, shape = 'line'),
      name = 'Ideal',
      hoverinfo = 'none'
    )
  plot <- plot %>%
    add_trace(
      r = self,
      theta = poles,
      name = paste("SSI Index:",round(self_index(wimp)$global[1],2)),
      marker = list(color = colors, size = 7, line = list(color = '#6F6BFF', width = 1.5)),
      fillcolor = 'rgba(204, 203, 248, 0.5)',
      line = list(width = 1, color = "#6F6BFF"),
      text = ~paste('<B>',construct,'</B>', '\nSelf:', round(self, 2), '\nIdeal:', round(ideal,2)),
      hoverinfo = 'text',
      hoverlabel=list(bgcolor = colors)
    )
  plot <- plot %>%
    layout(
      showlegend = FALSE,
      polar = list(
        radialaxis = list(
          visible = T,
          range = c(-1,1)
        )
      )
    )

  return(plot)

}

# SSI Heatmap -----------------------------------------------------------
#'
#' SSI Heatmap -- ssi_heatmap()
#'
#' @description A heat map representing the fit between I-actual and I-ideal as
#'              a function of the person's different cognitive states. It uses
#'              SSI Index in its calculations.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by  the \code{\link{importwimp}} function.
#'
#' @return A plotly heatmap.
#'
#' @author Alejandro Sanfeliciano
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#' ssi_heatmap(example.wimp)
#'

ssi_heatmap <- function(wimp){

  x <- .wimp_get_self(wimp)
  y <- .wimp_get_ideal(wimp)

  alpha.values <- seq(0, 1, by = 0.01)
  beta.values <- seq(0, 1, by = 0.01)

  sim_matrix <- outer(alpha.values, beta.values, Vectorize(function(alpha, beta) {
    .sim_index(x, y, alpha = alpha, beta = beta)
  }))

  plot <- plot_ly(
    x = alpha.values,
    y = beta.values,
    z = t(sim_matrix),
    type = "heatmap",
    colorscale = list(c(0, "#F52722"), c(0.5, "white"), c(1, "#A5D610")),
    zmin = 0,
    zmax = 1,
    hovertemplate = '<b>Alpha:</b> %{x}<br><b>Beta:</b> %{y}<br><b>Similarity:</b> %{z}<extra></extra>',
    colorbar = list(
      title = '<b>SSI</b>',
      tickfont = list(size = 16),
      ticklen = 10
    )
  ) %>%
    layout(
      title = "",
      xaxis = list(
        title = list(
          text = "<b>Discrepancy Salience (Alpha)</b>",
          font = list(size = 25)
        ),
        tickfont = list(size = 18)
        ),
      yaxis = list(
        title = list(
          text = "<b>Aspiration Salience (Beta)</b>",
          font = list(size = 25)
        ),
        tickfont = list(size = 18)),
      shapes = list(
        list(
          type = "rect",
          x0 = -0.005,
          x1 = 1.005,
          y0 = -0.005,
          y1 = 1.005,
          line = list(color = "black", width = 2)),
        list(
          type = "line",
          x0 = 0,
          y0 = 1,
          x1 = 1,
          y1 = 0,
          line = list(
            color = "black",
            width = 1,
            dash = "dot"
          )
        )
      ),
      annotations = list(
        list(
          x = 0.5,
          y = 0.5,
          xref = "x",
          yref = "y",
          text = "+",
          showarrow = FALSE,
          font = list(color = "black", size = 20)
        )
      )
    ) %>%
    style(
      hoverlabel = list(
        bgcolor = 'rgba(255, 255, 255, 0.8)',
        bordercolor = 'black',
        font = list(size = 12)
      )
    )

  return(plot)
}

#' Hypothetical Scenarios Plot  -- hypo_plot()
#'
#' @description This function creates a scatter plot to show the results of the
#'              \code{\link{self_index}} function.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param text.size Scalar that modifies the text size. Default is 1.
#' @param center Establishes the centre of the frame. Use "data" to set the data
#'        to be framed and "origin" to set the origin to be in the centre. the default
#'        is "data".
#' @param ... additional arguments are passed from \code{\link{self_index}}
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
#' hypo_plot (example.wimp)

hypo_plot <- function(wimp, text.size = 1, ...) {

  hypo.matrix <- .hypo.matrix(wimp)
  ncol <- ncol(hypo.matrix)
  hypo.names <- colnames(hypo.matrix)[-c(1,ncol)]

  self_index_data <- self_index(wimp, ...)


  congruence <- self_index_data$construct[[2]]

  construct.color <- ifelse(
    congruence == "Congruent", "#A5D610",
    ifelse(congruence == "Discrepant", "#F52722",
           ifelse(congruence == "Undefined", "yellow", "#000000"))
  )

  # Set up data.frame
  df <- self_index_data$construct[c(4,3)]
  df <- data.frame (df, construct.color, hypo.names)


  # Row and col names for data.frame
  names(df) <- c("ideal", "self", "color", "construct")
  rownames(df) <- hypo.names

  # Plotting
  y_ref <- self_index_data[[1]][[1]]

  fig <- plot_ly(
    data = df,
    x = ~self,
    y = ~ideal
  ) %>%
    add_annotations(
      data = df,
      x = ~self,
      y = ~ideal,
      text = ~construct,
      hoverinfo = 'text',
      font = list(size = 15 * text.size),
      showarrow = FALSE,


      xanchor = ~ifelse(self < 0.15, 'left', ifelse(self > 0.85, 'right', 'center')),
      xshift = ~ifelse(self < 0.15, 5, ifelse(self > 0.85, -5, 0)),
      yanchor = ~ifelse(ideal > 0.9, 'top', 'bottom'),
      yshift = ~ifelse(ideal > 0.9, -5, 5)

    ) %>%
    add_markers(
      data = df,
      x = ~self,
      y = ~ideal,
      marker = list(color = ~color, size = 7, line = list(color = 'black', width = 1)),
      text = ~paste('<b>', construct, '</b>', '\nIdeal Similarity:', ideal, '\nSelf Similarity:', self),
      hoverinfo = 'text'
    ) %>%
    layout(
      xaxis = list(
        title = "SELF SIMILARITY",
        range = c(0,1),
        gridwidth = 0.5,
        zeroline = TRUE,
        zerolinecolor = "black",
        zerolinewidth = 2
      ),
      yaxis = list(
        title = "IDEAL SIMILARITY",
        range = c(0,1),
        gridwidth = 0.5,
        zeroline = TRUE,
        zerolinecolor = "black",
        zerolinewidth = 2
      ),
      showlegend = FALSE,

      shapes = list(
        # Área verde (por encima)
        list(
          type = "rect",
          x0 = 0, x1 = 1,
          y0 = y_ref, y1 = 1,
          xref = "x", yref = "y",
          fillcolor = "rgba(0, 255, 0, 0.1)",
          line = list(width = 0)
        ),
        # Área roja (por debajo)
        list(
          type = "rect",
          x0 = 0, x1 = 1,
          y0 = 0, y1 = y_ref,
          xref = "x", yref = "y",
          fillcolor = "rgba(255, 0, 0, 0.1)",
          line = list(width = 0)
        ),
        # Línea horizontal negra
        list(
          type = "line",
          x0 = 0, x1 = 1,
          y0 = y_ref, y1 = y_ref,
          xref = "x", yref = "y",
          line = list(color = "#A5D610", dash = "dash", width = 2)
        )
      )
    )

  return(fig)

}
