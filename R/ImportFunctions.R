## IMPORT FUNCTIONS ##

# Import Weigthed ImpGrid -------------------------------------------------

#' Import Weighted Implication Grid -- importwimp()
#'
#' @description Function to transform the data of a WimpGrid contained in an
#' Excel file into an S3 object of class wimp.
#'
#' @param path Path to the excel file on your computer. The file suffix has to be .xlsx.
#' @param sheet Number of the Excel sheet that contains the WimpGrid data.
#'
#' @return A wimp S3 object.
#'
#' @export
#'
#' @import readxl


importwimp <- function(path, sheet = 1){

  wimp <- list()
  class(wimp) <- c("wimp","list")

  xlsx <- readxl::read_excel(path, sheet = sheet)

  n.constructs <- dim(xlsx)[1]

  # Scale -------------------------------------------------------------------
  scale.min <- as.numeric(names(xlsx)[1])
  scale.max <- as.numeric(names(xlsx)[n.constructs + 3])
  scale.center <- (scale.min + scale.max)/2
  scale <- c(scale.min,scale.max)
  names(scale) <- c("min","max")
  wimp$scale <- scale

  # Constructs --------------------------------------------------------------
  left.poles <- as.vector(xlsx[,1])[[1]]
  right.poles <- as.vector(xlsx[,n.constructs + 3])[[1]]
  constructs <- paste(left.poles,"-",right.poles,sep = " ")
  wimp$constructs[[1]] <- left.poles
  wimp$constructs[[2]] <- right.poles
  wimp$constructs[[3]] <- constructs
  names(wimp[["constructs"]]) <- c("left.poles","right.poles","constructs")

  # Self vector -------------------------------------------------------------
  xlsx.scores <- as.numeric(as.matrix(xlsx[,1:n.constructs+1]))
  direct.scores <- matrix(xlsx.scores,ncol=n.constructs,nrow=n.constructs)
  direct.self <- as.numeric(diag(direct.scores))
  standarized.self <- (direct.self - (scale.center * rep(1,n.constructs))) / (0.5 * (scale.max - scale.min))
  wimp$self[[1]] <- direct.self
  wimp$self[[2]] <- standarized.self
  names(wimp$self) <- c("direct","standarized")


  # Ideal vector ------------------------------------------------------------
  direct.ideal <- as.numeric(as.vector(xlsx[,n.constructs + 2])[[1]])
  standarized.ideal <- (direct.ideal - (scale.center * rep(1,n.constructs))) / (0.5 * (scale.max - scale.min))

  wimp$ideal[[1]] <- direct.ideal
  wimp$ideal[[2]] <- standarized.ideal
  names(wimp$ideal) <- c("direct","standarized")

  # Hypothetical vector -----------------------------------------------------
  standarized.hypothetical <- mapply(.calc.hypo, standarized.self,standarized.ideal)
  direct.hypothetical <- (scale.center * rep(1,n.constructs)) + (standarized.hypothetical * (0.5 * (scale.max - scale.min)))

  wimp$hypothetical[[1]] <- direct.hypothetical
  wimp$hypothetical[[2]] <- standarized.hypothetical
  names(wimp$hypothetical) <- c("direct","standarized")

  # Construct Categories------------------------------------------------------

  ind.dilemmatics <- which(standarized.ideal == 0)
  ind.undefined <- which(standarized.self == 0)
  ind.congruent <- which((standarized.ideal/standarized.self > 0) & !is.infinite(standarized.ideal/standarized.self))
  ind.discrepant <- which((standarized.ideal/standarized.self < 0) & !is.infinite(standarized.ideal/standarized.self))

  self.poles <- mapply(.self.poles, standarized.self,left.poles,right.poles)

  wimp$constructs$self.poles <- self.poles
  wimp$constructs$congruents <- ind.congruent
  wimp$constructs$discrepants <- ind.discrepant
  wimp$constructs$dilemmatics <- ind.dilemmatics
  wimp$constructs$undefined <- ind.undefined

  # Scores ------------------------------------------------------------------
  imp.matrix <- t((direct.scores - (scale.center * matrix(rep(1,n.constructs * n.constructs),ncol = n.constructs))) / (0.5 * (scale.max - scale.min)))

  num.weight.matrix <- imp.matrix - matrix(standarized.self,nrow = n.constructs,ncol = n.constructs,byrow = TRUE)
  den.weigth.matrix <- matrix(standarized.hypothetical,nrow = n.constructs,ncol = n.constructs) - matrix(standarized.self,nrow = n.constructs,ncol = n.constructs)
  weight.matrix <- num.weight.matrix / den.weigth.matrix

  wimp$scores[[1]] <- direct.scores
  wimp$scores[[2]] <- imp.matrix
  wimp$scores[[3]] <- weight.matrix
  names(wimp$scores) <- c("direct","implications","weights")

  # Function return ---------------------------------------------------------
  return(wimp)
}
