### HIDE FUNCTIONS ###


#' @importFrom dplyr select arrange case_when
#' @importFrom magrittr %>%
#' @importFrom igraph degree layout_with_graphopt layout_in_circle
#' @importFrom igraph layout_as_tree layout_with_mds layout_on_grid
#' @importFrom igraph betweenness closeness



# Threshold Function -----------------------------------------------
.thr <- function(x, method = "saturation"){

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
  # New format: values are already normalized; preserve as-is.
  # Optionally, you could implement pole-swapping based on ideal sign here.
  return(wimp)
}

# Hypothetical Situations Vector calculation ------------------------------------------------
.calc.hypo <- function(self, ideal) {
  # Handle NA values
  if(is.na(self) || is.na(ideal)) return(NA_real_)
  
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
  if(infer == "impact dynamics"){return("IMPACT")}
}

# Self Construct detection---------------------------
.self.poles <- function(self, l.pole, r.pole){

  construct <- paste(l.pole, "-", r.pole)

  # Handle NA or missing values
  if(is.na(self)) return(construct)
  
  # self is standardized value: negative=left pole, positive=right pole, 0=construct
  if(self < 0) return(l.pole)
  
  if(self > 0) return(r.pole)
  
  if(self == 0) return(construct)
  
  return(construct)
}

# Hypothetical Self matrix-----------------------------
.hypo.matrix <- function(wimp){

  imp.matrix <- wimp$scores$implications
  hypo.vector <- wimp$hypothetical$normalized
  self.vector <- wimp$self$normalized
  ideal.vector <- wimp$ideal$normalized

  constructs <- wimp$constructs$constructs
  left.poles <- wimp$constructs$left.poles
  right.poles <- wimp$constructs$right.poles
  hypo.poles <- mapply(.self.poles, hypo.vector,left.poles,right.poles)
  hypo.names <- hypo.poles

  hypo.matrix <- t(imp.matrix)
  diag(hypo.matrix) <- hypo.vector

  result <- cbind(self.vector,hypo.matrix,ideal.vector)

  colnames(result) <- c("SELF", hypo.names, "IDEAL")
  rownames(result) <- constructs

  return(result)
}

# Construct colors ------------------------------------------------------------
.construct.colors <- function(wimp, mode){
  col.sel <- .color.selection(mode)
  stopifnot(!is.null(wimp$vertices))
  n <- nrow(wimp$vertices)
  labels <- paste(wimp$vertices$lpole, "-", wimp$vertices$rpole)
  colors <- rep(NA_character_, n)
  self <- wimp$vertices$self
  ideal <- wimp$vertices$ideal
  idx_dilem <- ideal == 0
  idx_undef <- self == 0 & !idx_dilem
  idx_congr <- sign(self) == sign(ideal) & self != 0 & ideal != 0
  idx_discr <- !(idx_dilem | idx_undef | idx_congr)
  colors[idx_discr] <- col.sel[1]
  colors[idx_congr] <- col.sel[2]
  colors[idx_undef] <- col.sel[3]
  colors[idx_dilem] <- col.sel[4]
  res <- matrix(colors, ncol = 1)
  rownames(res) <- labels
  colnames(res) <- "color"
  res
}

# Merge two wimps --------------------------------------------------------------
.merge.wimp <- function(wimp1, wimp2){
  df1 <- data.frame(
    Construct = paste(wimp1$vertices$lpole, "-", wimp1$vertices$rpole),
    lpoles = wimp1$vertices$lpole,
    rpoles = wimp1$vertices$rpole,
    index1 = seq_len(nrow(wimp1$vertices))
  )
  df2 <- data.frame(
    Construct = paste(wimp2$vertices$lpole, "-", wimp2$vertices$rpole),
    lpoles = wimp2$vertices$lpole,
    rpoles = wimp2$vertices$rpole,
    index2 = seq_len(nrow(wimp2$vertices))
  )
  merge(df1, df2, by = 1:3, sort = FALSE)
}

# Compatibility merge wimps --------------------------------------------------------------
.compatibility.merge.wimp <- function(wimp1,wimp2){
  m <- nrow(.merge.wimp(wimp1,wimp2))
  n1 <- nrow(wimp1$vertices); n2 <- nrow(wimp2$vertices)
  if(m == 0) return("Incompatibility")
  if(m == n1 && m == n2) return("Full Compatibility")
  "Partial Compatibitily"
}

# Tversky Similarity function---------------------------------------------------
.sim_index <- function(x,y,alpha = .5, beta = .5){

  s.vec <- ifelse(y^2 > 0.25,
                  (-(x - y)^2) / (y^2) + 1,
                  (-(x - y)^2) / ((1 - abs(y)^2) + 1))

  vec.int <- s.vec[which(x*y > 0)]
  vec.x.minus.y <- x[which(x*y <= 0)]
  vec.y.minus.x <- y[which(x*y <= 0)]

  int.xy <- sum(abs(vec.int))
  x.minus.y <- sum(abs(vec.x.minus.y))
  y.minus.x <- sum(abs(vec.y.minus.x))

  sim.ratio <- int.xy / (int.xy + alpha * x.minus.y + beta * y.minus.x)
  return(sim.ratio)

}
# Isolate Construct ------------------------------------------------------------
.isolate.construct <- function(wimp, construct){
  # Keep only nodes connected to the specified construct (by row/col non-zero)
  stopifnot(!is.null(wimp$global$wmatrix))
  wmatrix <- wimp$global$wmatrix
  wmatrix[-construct, -construct] <- 0
  m <- which(apply(wmatrix, 1, function(row) any(row != 0)))
  n <- which(apply(wmatrix, 2, function(col) any(col != 0)))
  keep <- sort(unique(c(m, n)))
  # Reindex
  map <- match(keep, keep)
  new_vertices <- wimp$vertices[keep, , drop = FALSE]
  new_vertices$id <- seq_len(nrow(new_vertices))
  new_wmatrix <- wmatrix[keep, keep, drop = FALSE]
  # Build edges from filtered matrix
  edges <- NULL
  if(!is.null(new_wmatrix)){
    e_from <- c(); e_to <- c(); e_w <- c(); e_id <- c()
    for(i in seq_len(nrow(new_wmatrix))) for(j in seq_len(ncol(new_wmatrix))){
      w <- new_wmatrix[i, j]; if(i != j && !is.na(w) && w != 0){
        e_from <- c(e_from, i); e_to <- c(e_to, j); e_w <- c(e_w, w); e_id <- c(e_id, paste0(i, "t", j))
      }
    }
    edges <- data.frame(id = e_id, from = e_from, to = e_to, weight = e_w, stringsAsFactors = FALSE)
  }
  wimp$vertices <- new_vertices
  wimp$global$wmatrix <- new_wmatrix
  wimp$edges <- edges
  wimp
}

# Similarity function-----------------------------------------------------------
.sim <- function(s,i){

  result <- 1 - ((s - i)^2 / 4)
  return(result)
}

# Impact function---------------------------------------------------------------
.impact <- function(w,s,i){
  s.1 <- .thr(s + w)
  result <- .sim(s.1,i) - .sim(s,i)
  return(result)
}
