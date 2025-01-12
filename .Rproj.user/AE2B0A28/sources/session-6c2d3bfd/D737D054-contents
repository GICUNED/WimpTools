## System Dynamics Functions ##

# Scenario Matrix -----------------------------------------------------------

#' Scenario Matrix -- scenariomatrix()
#'
#' @description This function generates a scenario matrix based on the provided
#'              WimpGrid object and a vector of actions. It models potential
#'              outcomes iteratively based on user-defined parameters, allowing
#'              for scenario analysis and decision-making.
#'
#' @param wimp A "wimp" S3 object.
#' @param infer A character string indicating the type of inference to be used
#'              in the scenario analysis. Possible values include "self dynamics"
#'              or "impact dynamics". Default is "self dynamics"
#' @param thr A character string specifying the threshold method to be applied
#'            to the constructs. Supported options include "saturation","tanh"
#'            or "none". Default is "saturation"
#' @param act.vector A numeric vector representing the actions or interventions
#'                   to be applied to the constructs in the WimpGrid.
#' @param max.iter An integer specifying the maximum number of iterations for
#'                 the scenario analysis. Defaults is 5.
#' @param e A numeric value representing the convergence threshold. Smaller
#'          values lead to stricter convergence criteria. Default is 0.0001
#' @param stop.iter An integer specifying the number of consecutive iterations
#'                  with no significant change required to declare convergence.
#'                  Default is 2
#' @param exclude.dilemmatics A logical value. If TRUE, dilemmatic constructs
#'                            are excluded from the scenario matrix calculations.
#'                            Default is TRUE
#'
#' @return A "scn" S3 class object.
#' @export
#'
#' @examples
#'
#' # Activation vector (example)
#' act.vector <- c(0.1, -0.2, 0.3, -0.1, 0)
#'
#' # Generate a scenario matrix
#' scenariomatrix(wimp = example.wimp, act.vector = act.vector)
#'

scenariomatrix <- function(wimp, act.vector = NA, infer = "self dynamics",
                           thr = "saturation", max.iter = 5, e = 0.0001,
                           stop.iter = 2, exclude.dilemmatics = FALSE){


  if(!inherits(wimp,"wimp")){
    stop("The Weighted Implication Grid must be class wimp.")
  }
  if( ncol(wimp[[6]][[3]]) != length(act.vector) && infer != "impact dynamics"){
    stop("The weight matrix and the activation vector must have compatible dimensions.")
  }

  ideal <- wimp$ideal$standarized
  swap.vector <- ideal/abs(ideal)
  swap.vector[is.nan(swap.vector)] <- 1

  if(infer == "self dynamics"){
  act.vector <- act.vector * swap.vector
  }

  wimp <- .align.wimp(wimp,exclude.dilemmatics = exclude.dilemmatics)
  dim <- length(wimp$constructs$constructs)
  scene.matrix <- t(matrix(wimp[[3]][[2]]))
  trans.matrix <- t(wimp$scores$weights)
  next.matrix <- trans.matrix

  n <- 1
  i <- 0

  while(n <= max.iter && i <= stop.iter){

    if(infer == "self dynamics"){
      next.iter <- scene.matrix[n,] + t(act.vector)

      next.iter <- mapply(.thr, next.iter, thr)

      delta.iter <- next.iter - scene.matrix[n,]
      scene.matrix <- rbind(scene.matrix, next.iter)

      act.vector <- trans.matrix %*% delta.iter
    }

    if(infer == "impact dynamics"){
      if(n == 1){scene.matrix <- t(rep(0,dim))}
      if(exclude.dilemmatics == FALSE){
        n.matrix <- next.matrix
        n.matrix[.which.dilemmatics(wimp),] <- 0
      }else{
        n.matrix <- next.matrix
      }
      sum.columns <- t(n.matrix) %*% rep(1,nrow(trans.matrix))
      next.iter <- t(sum.columns)
      scene.matrix <- rbind(scene.matrix, next.iter)
      next.matrix <- trans.matrix %*% next.matrix
    }

    e.iter <- mean(abs(next.iter - scene.matrix[n,]))

    if(e.iter < e){
      i <- i + 1
    }else{
      i <- 0
    }
    n <- n + 1
  }

  rownames(scene.matrix) <- paste("iter", 0:(n-1))
  colnames(scene.matrix) <- wimp[[2]][[3]]

  if(n < max.iter){
    convergence <- n - (stop.iter + 1)
  }else{
    convergence <- NA
  }

  scene.list <- list()

  scene.list$values <- scene.matrix
  scene.list$convergence <- convergence
  scene.list$constructs <- wimp$constructs
  scene.list$self[[1]] <- wimp$self[[2]]
  scene.list$self[[2]] <- wimp$ideal[[2]]
  scene.list$weights <- wimp[[6]][[3]]
  scene.list$method$infer <- infer
  scene.list$method$threeshold <- thr

  names(scene.list$self) <- c("self","ideal")

  class(scene.list) <- c("scn","list")

  return(scene.list)
}

# PCSD -----------------------------------------------------------------

#' Personal Constructs System Dynamics plot -- pcsd()
#'
#' @description Interactive line plot of personal constructs system dinamics.
#' Show \code{\link{scenariomatrix}} values expressed in terms of distance to
#' Ideal-Self for each personal construct across the mathematical iterations.
#'
#' @param scn A "scn" S3 class object, the result of the \code{\link{scenariomatrix}} function.
#'            This object contains the scenario matrix and relevant metadata.
#' @param vline An optional numeric value indicating a specific iteration to highlight
#'              with a vertical line on the plot. Defaults is NA.
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
#'   wimp = example.wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act.vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Plot the dynamics
#' pcsd(example_scn)
#'


pcsd <- function(scn, vline = NA){

  poles <- scn$constructs$constructs
  dim <- length(poles)
  infer <- scn$method$infer
  iter <- nrow(scn$values)


  self.vector <- scn$self[[1]]
  self.matrix <- matrix(self.vector, ncol = length(self.vector),
                        nrow = iter, byrow = TRUE)

  res <- scn$values


  x <- c(0:(iter -1))
  y <- c(0:length(poles))
  y <- as.character(y)

  if(infer == "self dynamics"){
    df <- data.frame(x, (res - self.matrix))
  }
  if(infer == "impact dynamics"){
    df <- data.frame(x, (res/dim))
  }

  max.value.df <- max(abs(df[,-1])) + 0.05 * max(abs(df[,-1]))

  colnames(df) <- y

  fig <- plotly::plot_ly(df, x = ~x, y = df[,2], name = poles[1],
                         type = 'scatter',
                         mode = 'lines+markers',line = list(shape = "spline"))  # Build PCSD with plotly.

  for (n in 3:(length(poles)+1)) {
    fig <- fig %>% plotly::add_trace(y = df[,n], name = poles[n-1],
                                     mode = 'lines+markers',
                                     line = list(shape = "spline"))
  }
  fig <- fig %>% plotly::layout(
    xaxis = list(
      title = "ITERATIONS"
    ),
    yaxis = list(
      title = .label.y(infer),
      range = c(-max.value.df,max.value.df)
    )
  )
  fig <- fig %>% plotly::layout(legend=list(
    title=list(text='<b>PERSONAL CONSTRUCTS</b>')
  )
  )

  fig <- fig %>% add_lines(
    x = vline,
    y = c(-max.value.df,max.value.df),
    line = list(
      color = "grey",
      dash = "dot"
    ),
    inherit = FALSE,
    showlegend = FALSE
  )

  fig
}

# AUC Index ---------------------------------------------------------------

#' PCSD AUC Index -- auc_index()
#'
#' @description This function calculates the area under the PCSD curve for each
#' construct.
#'
#' @param scn A "scn" S3 class object, the result of the \code{\link{scenariomatrix}} function.
#'            This object contains the scenario matrix and relevant metadata.
#' @param std A logical value. If \code{TRUE}, standardizes the AUC results.
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
#'   wimp = example.wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act.vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate AUC index for each construct
#' auc_index(example_scn,std = TRUE)
#'

auc_index <- function(scn, std = TRUE){

  poles <- scn$constructs$constructs
  iter <- nrow(scn$values)
  self.vector <- scn$self$self
  self.matrix <- matrix(self.vector, ncol = length(self.vector),
                        nrow = iter, byrow = TRUE)
  std.coef <- ifelse(std, iter, 1)

  res <- scn$values

   if(scn$method$infer == "self dynamics"){
     res <- res - self.matrix
    }

  matrix <- matrix(ncol= length(poles), nrow = 1)

  for (n in 1:length(poles)) {
    matrix[,n] <- MESS::auc(c(0:(iter - 1)), res[,n], type = "spline")/std.coef
  }

  result <- t(matrix)

  rownames(result) <- poles
  colnames(result) <- "AUC"

  return(result)
}

# PCSD Stability Index ----------------------------------------------------

#' PCSD Stability Index -- stability_index()
#'
#' @description This function returns the standard deviation for each
#' construct over the mathematical iterations of the PCSD.
#'
#' @param scn A "scn" S3 class object, the result of the \code{\link{scenariomatrix}} function.
#'            This object contains the scenario matrix and relevant metadata.
#'
#' @return Returns a vector with the standard deviation of each of the
#' constructs.
#'
#' @importFrom stats sd
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example.wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act.vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate stability index for each construct
#' stability_index(example_scn)
#'

stability_index <- function(scn){

  poles <- scn$constructs$constructs
  iter <- nrow(scn$values)
  self.vector <- scn$self$self
  self.matrix <- matrix(self.vector, ncol = length(self.vector),
                        nrow = iter, byrow = TRUE)

  res <- scn$values

  result <- apply(res, 2, sd)                                                   # Calculate SD for each construct.

  result <- matrix(result)
  rownames(result) <- poles                                                     # Name vector's elements.
  colnames(result) <- "Standard Deviation"

  return(result)
}

# PCSD Summary ------------------------------------------------------------

#' PCSD summary -- pcsd_summary()
#'
#' @description This function returns a summary of the PCSD. It informs us the
#' initial and final value of each construct and the difference between them.
#'
#' @param scn A "scn" S3 class object, the result of the \code{\link{scenariomatrix}} function.
#'            This object contains the scenario matrix and relevant metadata.
#'
#' @return Returns a matrix with the PCSD summary.
#'
#'
#' @export
#'
#' @examples
#' # Example Scenario Matrix (scn)
#' example_scn <- scenariomatrix(
#'   wimp = example.wimp,
#'   infer = "self dynamics",
#'   thr = "saturation",
#'   act.vector = c(0.1, -0.2, 0.3, -0.1, 0)
#' )
#'
#' # Calculate PCSD summary for each construct
#' pcsd_summary(example_scn)
#'

pcsd_summary <- function(scn){


  poles <- scn$constructs$constructs
  iter <- nrow(scn$values)
  self.vector <- scn$self$self
  self.matrix <- matrix(self.vector, ncol = length(self.vector),
                        nrow = iter, byrow = TRUE)

  res <- scn$values
  average <- colMeans(res)

  result <- res[c(1,iter),]                                                     # Extract the first vector and the last vector from the iteration matrix.
  result <- t(result)
  result <- cbind(result, average, result[,2] - result[,1])                              # Calculate the difference between the first vector and the last one and add it to the results.

  rownames(result) <- poles
  colnames(result) <- c("Initial", "Final", "Average", "Difference")

  return(result)
}
