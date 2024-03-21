### HIDE FUNCTIONS ###

# Lineal Threshold Function -----------------------------------------------

.thr <- function(x, method){

  if(method == "none"){
    result <- x
  }

  if(method == "saturation"){
    if(x <= -1){ result <- -1}
    if(-1 < x && x < 1){ result <- x}
    if(x >= 1){ result <- 1}
  }

  if(method == "tanh"){
    result <- tanh(x)
  }

  return(result)
}


# Color function ----------------------------------------------------------

.color.selection <- function(x){                                                # Order: c(discrepant, congruent, undefined , dilemmatic)

  if(x == "red/green"){
    res <- c("#F52722","#A5D610","grey","yellow")
  }

  if(x == "grey scale"){
    res <- c("#808080","#ffffff","#f2f2f2","#e5e5e5")
  }

  return(res)
}


# Align wimp function -----------------------------------------------------

.align.wimp <- function(wimp, exclude.dilemmatics = TRUE){

  ideal <- wimp$ideal[[2]]

  swap.indeces <- which(ideal < 0)
  dil.indeces <- which(ideal == 0)

# Scale transformation

  wimp$scale[1] <- -1
  wimp$scale[2] <- 1

# Constructs transformation

  old.left.poles <- wimp$constructs$left.poles
  old.right.poles <- wimp$constructs$right.poles

  left.poles <- old.left.poles
  left.poles[swap.indeces] <- old.right.poles[swap.indeces]

  right.poles <- old.right.poles
  right.poles[swap.indeces] <- old.left.poles[swap.indeces]

  if(exclude.dilemmatics && length(dil.indeces) > 0){
    left.poles <- left.poles[-dil.indeces]
    right.poles <- right.poles[-dil.indeces]
  }

  constructs <- paste(left.poles,"—",right.poles,sep = "")

  wimp$constructs$left.poles <- left.poles
  wimp$constructs$right.poles <- right.poles
  wimp$constructs$constructs <- constructs

# Self transformation

  self <- wimp$self$standarized
  self[swap.indeces] <- self[swap.indeces] * -1

  if(exclude.dilemmatics && length(dil.indeces) > 0){
    self <- self[-dil.indeces]
  }

  wimp$self$direct <- self
  wimp$self$standarized <- self

# Ideal transformation

  ideal <- wimp$ideal$standarized
  ideal[swap.indeces] <- ideal[swap.indeces] * -1

  if(exclude.dilemmatics && length(dil.indeces) > 0){
    ideal <- ideal[-dil.indeces]
  }

  wimp$ideal$direct <- ideal
  wimp$ideal$standarized <- ideal

# Hypothetical transformation

  hypothetical <- wimp$hypothetical$standarized
  hypothetical[swap.indeces] <- hypothetical[swap.indeces] * -1

  if(exclude.dilemmatics && length(dil.indeces) > 0){
    hypothetical <- hypothetical[-dil.indeces]
  }

  wimp$hypothetical$direct <- hypothetical
  wimp$hypothetical$standarized <- hypothetical

# Scores transformation

  implications <- wimp$scores$implications
  weights <- wimp$scores$weights

  implications[swap.indeces,] <- implications[swap.indeces,] * -1
  implications[,swap.indeces] <- implications[,swap.indeces] * -1

  weights[swap.indeces,] <- weights[swap.indeces,] * -1
  weights[,swap.indeces] <- weights[,swap.indeces] * -1

  if(exclude.dilemmatics && length(dil.indeces) > 0){
    implications <- implications[-dil.indeces,]
    implications <- implications[,-dil.indeces]

    weights <- weights[-dil.indeces,]
    weights <- weights[,-dil.indeces]
  }

  wimp$scores$direct <- implications
  wimp$scores$implications <- implications
  wimp$scores$weights <- weights

# OpenRepGrid transformation

# Return
  return(wimp)
}

# Hypothetical calculation ------------------------------------------------

.calc.hypo <- function(self, ideal) {
  if (self != 0) {
    return(self / (-1 * abs(self)))

  } else if (self == 0 && !(0 %in% ideal)) {
    return(ideal / abs(ideal))

  } else if (self == 0 && (0 %in% ideal)) {
    return(1)
  }
}

# Dilemmatics detection ---------------------------------------------------

.which.dilemmatics <- function(wimp){
  ideal <- wimp$ideal[[2]]
  dil.indeces <- which(ideal == 0)
  return(dil.indeces)
}


# PCSD Y-Axis label -------------------------------------------------------

.label.y <- function(infer){
  if(infer == "self dynamics"){return("SELF DIFFERENTIAL")}
  if(infer == "adjustment dynamics"){return("ADJUSTMENT IMPACT")}
}

# Self Construct detection

.self.poles <- function(self,l.pole,r.pole){
  construct <- paste(l.pole,"—",r.pole)
  if(self > 0){return(r.pole)}
  if(self < 0){return(l.pole)}
  if(self == 0){return(construct)}
}
