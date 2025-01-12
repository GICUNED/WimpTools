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

n.congruents <- length(wimp$constructs$congruents)
n.discrepants <- length(wimp$constructs$discrepants)
n.dilemmatics <- length(wimp$constructs$dilemmatics)
n.undefined <- length(wimp$constructs$undefined)

n <- length(wimp$self$standarized)

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
#'        are allowed and "sim" for Tversky similarity. Default is "sim".
#' @param rc Use Cohen's rc which is invariant to construct reflection. Default is TRUE.
#' @param alpha Alpha value for Tversky similarity calculation. Default is .5.
#' @param beta Beta value for Tversky similarity calculation. Default is .5.
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

self_index <- function(wimp, method = "sim", rc = TRUE, alpha = .5, beta = .5){

  result <- list()

  congruence <- rep("Undefined",length(wimp$self$standarized))
  congruence[wimp$constructs$congruents] <- "Discrepant"
  congruence[wimp$constructs$undefined] <- "Congruent"
  congruence[wimp$constructs$discrepants] <- "Congruent"
  congruence[wimp$constructs$dilemmatics] <- "Undefined"

  hypo.matrix <- .hypo.matrix(wimp)
  ncol <- ncol(hypo.matrix)
  hypo.names <- colnames(hypo.matrix)[-c(1,ncol)]
  hypo.names <- paste("Totally", hypo.names, sep = " ")

  rc.text <- "no rc"
  if(rc){
    hypo.matrix <- rbind(hypo.matrix, -hypo.matrix)
    rc.text <- "rc"
    }

  if(method == "sim"){

    self.vector <- hypo.matrix[,1]
    ideal.vector <- hypo.matrix[,ncol]

    discrepants <- wimp$constructs$discrepants
    congruents <- wimp$constructs$congruents

    self.cor <- apply(hypo.matrix[,-c(1, ncol(hypo.matrix))], 2, function(col) .sim_index(self.vector, col, alpha = alpha, beta = beta))
    ideal.cor <- apply(hypo.matrix[,-c(1, ncol(hypo.matrix))], 2, function(col) .sim_index(ideal.vector, col, alpha = alpha, beta = beta))

    self.ideal.cor <- .sim_index(self.vector,ideal.vector, alpha = alpha, beta = beta)
    self.hypo.cor <- mean(self.cor)
    ideal.hypo.cor <- mean(ideal.cor)

    ideal.hypo.congruents.cor <- mean(ideal.cor[congruents])
    ideal.hypo.discrepants.cor <- mean(ideal.cor[discrepants])

  }
  if(!( method == "sim" | method == "cos" )){
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

  df.construct <- data.frame(cbind(hypo.names,congruence,round(as.numeric(self.cor),4),round(as.numeric(ideal.cor),4)))
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

  self <- wimp$self$standarized
  self <- c(self,self[1])

  ideal <- wimp$ideal$standarized
  ideal <- c(ideal,ideal[1])

  r.poles <- wimp$constructs$right.poles
  l.poles <- wimp$constructs$left.poles
  poles <- paste(r.poles," (",l.poles,")", sep="")
  poles <- c(poles,poles[1])

  construct<- wimp$constructs$constructs
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
      marker = list(color = colors, size = 7, line = list(color = '#FA9D13', width = 1.5)),
      fillcolor = 'rgba(255, 217, 125, 0.5)',
      line = list(width = 1, color = "#FA9D13"),
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

# Adjustment Heatmap -----------------------------------------------------------
#'
#' Adjustment Heatmap -- ssi_heatmap()
#'
#' @description A heat map representing the fit between I-actual and I-ideal as
#'              a function of the person's different cognitive states. It uses
#'              Tversky's concept of similarity in its calculations.
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

  x <- wimp$self$standarized
  y <- wimp$ideal$standarized

  alpha.values <- seq(0, 1, by = 0.01)
  beta.values <- seq(0, 1, by = 0.01)

  sim_matrix <- outer(alpha.values, beta.values, Vectorize(function(alpha, beta) {
    .sim_index(x, y, alpha = alpha, beta = beta)
  }))

  plot <- plot_ly(
    x = alpha.values,
    y = beta.values,
    z = sim_matrix,
    type = "heatmap",
    colorscale = list(c(0, "#F52722"), c(0.5, "white"), c(1, "#A5D610")),
    zmin = 0,
    zmax = 1,
    hovertemplate = '<b>Alpha:</b> %{x}<br><b>Beta:</b> %{y}<br><b>Adjustment:</b> %{z}<extra></extra>'
  ) %>%
    layout(
      title = "",
      xaxis = list(title = "Attention to Self Discrepancies (Alpha)"),
      yaxis = list(title = "Attention to the Desired Change (Beta) "),
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
