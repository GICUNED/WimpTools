## DIGRAPH FUNCTIONS ##


# Self Digraph ------------------------------------------------------------

#' Selfdigraph -- digraph()
#'
#' @description A digraph that represents the self of the person being assessed
#'              on the basis of its constructs and the relationships between them.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#'
#' @param vertex.vector Vector defining the value of each of the vertices of the digraph.
#'        Default is the value of the standarized self from the wimp object.
#' @param ideal.vector Vector defining the ideal value of each of the vertices of the digraph.
#'        Default is the value of the standarized ideal self from the wimp object.
#' @param width digraph width.
#' @param height Digraph heigth.
#' @param color Color palette to be used. The options are "red/green" and "grey scale". Default is "red/green".
#' @param layout Layout with which the digraph will be displayed. The options
#'        are: "circle", "rtcircle", "tree", "graphopt", "mds" and "grid". Default is "graphopt".
#' @param show Logical vector defining which constructs to display. By default all constructs are shown.
#' @param hide.direct If TRUE, hide direct relationship between nodes of the graph. Default is FALSE.
#' @param areas Logical; if TRUE, draw coloured areas that group nodes by category. Default is TRUE.
#' @param pad_side Numeric; padding in pixels used for the areas drawn around the nodes. Default is 36.
#' @param rounding Numeric; radius used to round the area corners (0 for sharp corners). Default is 0.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A digraph made with visNetwork
#'
#' @import visNetwork
#' @importFrom htmlwidgets JS
#' @importFrom jsonlite toJSON
#' @export
#'
#' @examples
#'
#'  digraph(su.wimp)
#'

digraph <- function(wimp, vertex.vector = NA, ideal.vector = NA, width="100%",
                    height="1000px", color = "red/green", layout ="graphopt",
                    show = TRUE, hide.direct = FALSE,
                    areas = TRUE, pad_side = 36, rounding = 0){

  if(inherits(wimp,"wimp")){
    lpoles <- wimp$constructs[[1]]
    rpoles <- wimp$constructs[[2]]
    poles <- paste(lpoles, "-", rpoles)
    wmatrix <- wimp[["scores"]][["weights"]]
    self <- wimp$self[[2]]
    ideal <- wimp$ideal[[2]]
    cat_vec <- wimp$constructs$category}

  if(inherits(wimp,"scn")){
    lpoles <- wimp$constructs[[1]]
    rpoles <- wimp$constructs[[2]]
    poles <- paste(lpoles, "-", rpoles)
    wmatrix <- wimp$weights
    self <- wimp[["self"]][[1]]
    ideal <- wimp[["self"]][[2]]
    cat_vec <- wimp$constructs$category}

  if(!exists("cat_vec")){
    cat_vec <- rep(NA_character_, length(lpoles))
  }

  if(is.na(vertex.vector[1])){
    vertex.vector <-  self
  }else{
    vertex.vector <- vertex.vector
  }
    vertex.vector <- sapply(vertex.vector,.thr)

  if(is.na(ideal.vector[1])){
    ideal.vector <- ideal
  }else{
    ideal.vector <- ideal.vector
  }

  vertex.name <- c()
  n <- 1
  for (x in vertex.vector) {
    if(x < 0){vertex.name[n] <- lpoles[n] }
    else{
      if(x > 0){vertex.name[n] <- rpoles[n] }
      else{
        if(x == 0){vertex.name[n] <- poles[n]}
      }
    }
    n <- n + 1
  }

  congruency.vector <- vertex.vector/ideal.vector

  discrepant.color <- .color.selection(color)[1]
  congruent.color <- .color.selection(color)[2]
  undefined.color <- .color.selection(color)[3]
  dilemmatic.color <- .color.selection(color)[4]

  vertex.color <- sapply(congruency.vector, function(x) case_when(
    x < 0 && x != -Inf ~ discrepant.color,
    x > 0 && x != Inf ~ congruent.color,
    x == 0 ~ undefined.color,
    is.infinite(x) ~ dilemmatic.color,
    .default = dilemmatic.color)
  )

  vertex.group <- sapply(congruency.vector, function(x) case_when(
    x < 0 && x != -Inf ~ "Discrepant",
    x > 0 && x != Inf ~ "Congruent",
    x == 0 ~ "Undefined Self",
    is.infinite(x) ~ "Dilemmatic",
    .default = "Dilemmatic")
  )

  level.vector <- vertex.vector
  level.vector[level.vector == 0] <- 1
  level.vector <- level.vector * ideal.vector

  vertex.vadjust <- sapply(abs(vertex.vector), function(x) -4*x^2-25*x-42)

  vertex <- data.frame(
    id = 1:length(vertex.vector),
    label = vertex.name,
    group = vertex.group,
    category = as.character(cat_vec),
    size = 30 * abs(vertex.vector) + 20,
    shape = "dot",
    title = paste("<p><b>", poles,"</b><br>Self:",
                  round(vertex.vector,2),"<br>Ideal:",
                  round(ideal,2),"</p>"),
    color = vertex.color,
    shadow = TRUE,
    hidden = !show,
    font.size = 20,
    font.strokeWidth = 3,
    font.vadjust = vertex.vadjust
  )

  n <- 1

  for (i in vertex.vector){
    if(i != 0){
      direction.value <- i / abs(i)
      wmatrix[,n] <- wmatrix[,n] * direction.value
      wmatrix[n,] <- wmatrix[n,] * direction.value
    }
    n <- n + 1
  }
  n.vertex <- nrow(wmatrix)

  if(hide.direct){
    logical.dilemmatic <- ideal == 0

    wmatrix[wmatrix > 0] <- 0
    wmatrix[logical.dilemmatic,] <- 0
    wmatrix[,logical.dilemmatic] <- 0
  }
  weight.vector <- c()
  from.vector <- c()
  to.vector <- c()

  for (i in 1:n.vertex) {
    for (j in 1:n.vertex) {
      weight <- wmatrix[i, j]
      if (weight != 0) {
        weight.vector <- c(weight.vector,weight)
        from.vector <- c(from.vector,i)
        to.vector <- c(to.vector,j)
      }}}

  if(color!= "grey scale" ){
    edges.color <- sapply(weight.vector, function(x)
      ifelse(x > 0 , "grey","#CD5C5C"))
    edges.dashes <- FALSE
  }else{
    edges.color <- "grey"
    edges.dashes <- sapply(weight.vector, function(x)
      ifelse(x > 0 , FALSE,TRUE))
  }

  edge.curved <- logical(length(weight.vector))
  n <- 1
  for (N in 1:dim(wmatrix)[1]) {
    for (M in 1:dim(wmatrix)[1]) {
      if(wmatrix[M,N] != 0 && wmatrix[N,M] != 0){
        edge.curved[n] <- TRUE
      }
      if(wmatrix[N,M] != 0){
        n <- n + 1
      }}}

  edges <- data.frame(
    from = from.vector,
    to = to.vector,
    width = 2 * abs(weight.vector),
    arrows = "to",
    dashes = edges.dashes,
    smooth = edge.curved,
    color = edges.color,
    title = round(weight.vector, 2)
  )

  if(layout == "graphopt") {layout <- "layout_with_graphopt"}
  if(layout == "circle") {layout <- "layout_in_circle"}
  if(layout == "tree") {layout <- "layout_as_tree"}
  if(layout == "mds") {layout <- "layout_with_mds"}
  if(layout == "grid") {layout <- "layout_on_grid"}

  g <- if(layout == "rtcircle") {
    visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = "layout_as_tree", circular = TRUE) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 0, labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
  } else {
    visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = layout, randomSeed = 33) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 0, labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
  }

  if("category" %in% names(vertex) && !is.null(g$x$nodes)){
    if(NROW(g$x$nodes) == nrow(vertex)){
      g$x$nodes$category <- vertex$category
    }
  }

  if(areas){
    uniq_cats <- unique(stats::na.omit(vertex$category[vertex$category != ""]))
    if(length(uniq_cats) > 0){
      base_cols <- c("#6EA8FE","#72D6A0","#F7B267","#D985B9","#8FD3FE",
                     "#B5E48C","#FFD166","#CDB4DB","#A0C4FF","#FFAFCC")
      fill_vec   <- rep(base_cols, length.out = length(uniq_cats))
      stroke_vec <- fill_vec
      rgb_to_rgba <- function(hex, a){
        hex <- gsub("#","",hex)
        r <- strtoi(substr(hex,1,2),16)
        g <- strtoi(substr(hex,3,4),16)
        b <- strtoi(substr(hex,5,6),16)
        sprintf("rgba(%d,%d,%d,%.2f)", r,g,b,a)
      }
      fills   <- vapply(fill_vec,   rgb_to_rgba, character(1), a = 0.15)
      strokes <- vapply(stroke_vec, rgb_to_rgba, character(1), a = 0.65)

      js_hulls <- paste0(
        "function(ctx){\n",
        "  var net=this;\n",
        "  var PAD=", pad_side, ";\n",
        "  var ROUND=", rounding, ";\n",
        "  var LINE_W=2;\n\n",
        "  var CAT = ", jsonlite::toJSON(uniq_cats, auto_unbox=TRUE), ";\n",
        "  var FILL= ", jsonlite::toJSON(as.character(fills), auto_unbox=TRUE), ";\n",
        "  var STROK=", jsonlite::toJSON(as.character(strokes), auto_unbox=TRUE), ";\n\n",
        "  function isCCW(pts){var s=0;for(var i=0;i<pts.length;i++){var a=pts[i],b=pts[(i+1)%pts.length];s+=a.x*b.y-a.y*b.x;}return s>0;}\n",
        "  function lineIntersect(p1,d1,p2,d2){var det=d1.x*d2.y-d1.y*d2.x;if(Math.abs(det)<1e-9) return null;var t=((p2.x-p1.x)*d2.y-(p2.y-p1.y)*d2.x)/det;return {x:p1.x+t*d1.x,y:p1.y+t*d1.y};}\n",
        "  function convexHull(points){\n",
        "    if(points.length<=1) return points.slice();\n",
        "    var pts=points.slice().sort(function(a,b){return a.x!==b.x?a.x-b.x:a.y-b.y;});\n",
        "    function cross(o,a,b){return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);}\n",
        "    var lower=[],upper=[];\n",
        "    for(var i=0;i<pts.length;i++){while(lower.length>=2 && cross(lower[lower.length-2], lower[lower.length-1], pts[i])<=0) lower.pop(); lower.push(pts[i]);}\n",
        "    for(var j=pts.length-1;j>=0;j--){while(upper.length>=2 && cross(upper[upper.length-2], upper[upper.length-1], pts[j])<=0) upper.pop(); upper.push(pts[j]);}\n",
        "    upper.pop(); lower.pop(); return lower.concat(upper);\n",
        "  }\n",
        "  function offsetConvex(pts, d){\n",
        "    if(pts.length<3) return null;\n",
        "    var H=pts.slice(); if(!isCCW(H)) H.reverse();\n",
        "    var N=H.length, out=new Array(N);\n",
        "    for(var i=0;i<N;i++){\n",
        "      var p0=H[(i-1+N)%N], p1=H[i], p2=H[(i+1)%N];\n",
        "      var e0={x:p1.x-p0.x,y:p1.y-p0.y}, e1={x:p2.x-p1.x,y:p2.y-p1.y};\n",
        "      var l0=Math.hypot(e0.x,e0.y)||1, l1=Math.hypot(e1.x,e1.y)||1; e0.x/=l0; e0.y/=l0; e1.x/=l1; e1.y/=l1;\n",
        "      var n0={x:e0.y,y:-e0.x}, n1={x:e1.y,y:-e1.x};\n",
        "      var pB={x:p1.x+n0.x*d,y:p1.y+n0.y*d}, pC={x:p1.x+n1.x*d,y:p1.y+n1.y*d};\n",
        "      var q=lineIntersect(pB,e0,pC,e1);\n",
        "      out[i]= q ? q : {x:p1.x+(n0.x+n1.x)*d, y:p1.y+(n0.y+n1.y)*d};\n",
        "    }\n",
        "    return out;\n",
        "  }\n",
        "  function draw(ctx, pts, fill, stroke){\n",
        "    if(!pts || pts.length<3) return;\n",
        "    ctx.save();\n",
        "    ctx.globalCompositeOperation='source-over';\n",
        "    ctx.fillStyle=fill; ctx.strokeStyle=stroke; ctx.lineWidth=LINE_W;\n",
        "    if(ROUND<=0){\n",
        "      ctx.beginPath(); ctx.moveTo(pts[0].x,pts[0].y);\n",
        "      for(var i=1;i<pts.length;i++) ctx.lineTo(pts[i].x,pts[i].y);\n",
        "      ctx.closePath(); ctx.fill(); ctx.stroke();\n",
        "    } else {\n",
        "      ctx.lineJoin='round'; ctx.lineCap='round'; ctx.miterLimit=2;\n",
        "      ctx.beginPath(); var n=pts.length;\n",
        "      for(var i=0;i<n;i++){\n",
        "        var p0=pts[(i-1+n)%n], p1=pts[i], p2=pts[(i+1)%n];\n",
        "        var v1x=p1.x-p0.x,v1y=p1.y-p0.y,v2x=p2.x-p1.x,v2y=p2.y-p1.y;\n",
        "        var l1=Math.hypot(v1x,v1y)||1,l2=Math.hypot(v2x,v2y)||1;\n",
        "        var rr=Math.min(ROUND,0.45*l1,0.45*l2);\n",
        "        v1x/=l1; v1y/=l1; v2x/=l2; v2y/=l2;\n",
        "        var p1_in={x:p1.x-v1x*rr,y:p1.y-v1y*rr}, p1_out={x:p1.x+v2x*rr,y:p1.y+v2y*rr};\n",
        "        if(i===0) ctx.moveTo(p1_in.x,p1_in.y); else ctx.lineTo(p1_in.x,p1_in.y);\n",
        "        ctx.arcTo(p1.x,p1.y,p1_out.x,p1_out.y,rr);\n",
        "      }\n",
        "      ctx.closePath(); ctx.fill(); ctx.stroke();\n",
        "    }\n",
        "    ctx.restore();\n",
        "  }\n\n",
        "  var nodes=net.body.data.nodes.get();\n",
        "  var bycat={};\n",
        "  nodes.forEach(function(n){\n",
        "    if(n.hidden) return;\n",
        "    if(n.category==null) return;\n",
        "    var c=String(n.category);\n",
        "    if(c==='' || c==='NA') return;\n",
        "    (bycat[c]||(bycat[c]=[])).push(n.id);\n",
        "  });\n\n",
        "  for(var i=0;i<CAT.length;i++){\n",
        "    var c = CAT[i];\n",
        "    var ids = bycat[c];\n",
        "    if(!ids || ids.length===0) continue;\n",
        "    var pos = net.getPositions(ids);\n",
        "    var pts = ids.map(function(id){ return {x:pos[id].x, y:pos[id].y}; });\n\n",
        "    if(pts.length<=2){\n",
        "      var xs=pts.map(function(p){return p.x;}), ys=pts.map(function(p){return p.y;});\n",
        "      var minX=Math.min.apply(null,xs)-PAD, maxX=Math.max.apply(null,xs)+PAD;\n",
        "      var minY=Math.min.apply(null,ys)-PAD, maxY=Math.max.apply(null,ys)+PAD;\n",
        "      pts=[{x:minX,y:minY},{x:maxX,y:minY},{x:maxX,y:maxY},{x:minX,y:maxY}];\n",
        "      draw(ctx, pts, FILL[i], STROK[i]);\n",
        "    }else{\n",
        "      var hull=convexHull(pts);\n",
        "      var off = offsetConvex(hull, PAD);\n",
        "      draw(ctx, off, FILL[i], STROK[i]);\n",
        "    }\n",
        "  }\n",
        "}"
      )

      g <- g %>% visEvents(beforeDrawing = htmlwidgets::JS(js_hulls))
    }
  }

  g
}
# Ideal Digraph -----------------------------------------------------------

#' Ideal digraph -- idealdigraph()
#'
#' @description A digraph that represents the ideal self of the person being assessed
#'              on the basis of its constructs and the relationships between them.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param inc If TRUE, hide direct relationship between nodes of the graph. Default is FALSE.
#' @param ... additional arguments are passed from \code{\link{digraph}}
#'        function.
#' @param layout Layout with which the digraph will be displayed. The options
#'        are: "circle", "rtcircle", "tree", "graphopt", "mds" and "grid". Default is "circle".
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A digraph made with visNetwork
#'
#' @export
#'
#' @examples
#'
#' idealdigraph(example.wimp)
#'

idealdigraph <- function(wimp, inc=FALSE, layout = "circle", ...){

  ideal.vector <- wimp$ideal[[2]]
  plot <- digraph(wimp = wimp, hide.direct = inc, vertex.vector = ideal.vector, layout = layout, ...)
  return(plot)
}

# Simulation Digraph ---------------------------------------------------------

#' Simulation digraph -- simdigraph()
#'
#' @description A digraph that represents the hypothetical self of the person being assessed
#'              on the basis of its constructs and the relationships between them.
#'
#' @param scn A scenario matrix. It must be a "scn" S3 object
#'         from the \code{\link{scenariomatrix}} function.
#' @param niter Iteration displayed in the digraph. Default is 0.
#' @param ... Additional arguments are passed from \code{\link{digraph}}
#'        function.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A digraph made with visNetwork
#'
#' @export
#'
#' @examples
#'
#'  scn <- scenariomatrix(example.wimp, rep(1,5))
#'  simdigraph(scn, niter = 2)
#'

simdigraph <- function(scn, niter = 0, ...){

  vertex.vector <- scn[[1]][niter+1,]
  digraph(wimp = scn, vertex.vector = vertex.vector, ...)
}

# In - Out Digraph -------------------------------------------------------------

#' In - Out digraph -- inout_digraph()
#'
#' @description A digraph showing the in and out vertices for a given construct.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param iso.con Objective construct to analyse the in and out. Must be an
#'        integer representing its position in the list of constructs.
#' @param ... additional arguments are passed from \code{\link{digraph}}
#'        function.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A digraph made with visNetwork
#'
#' @export
#'
#' @examples
#'
#' idealdigraph(example.wimp)

inout_digraph <- function(wimp,iso.con, ...){

  wimp <- .isolate.construct(wimp,iso.con)

  wmatrix <- abs(wimp$scores$weights)
  center <- which.max( colSums(wmatrix) + rowSums(wmatrix))

  out.vertex <- which(colSums(wmatrix)!=0)
  in.vertex <- which(rowSums(wmatrix)!=0)
  inout.v <- intersect(in.vertex,out.vertex)

  out.v <- setdiff(out.vertex,in.vertex)
  in.v <- setdiff(in.vertex,out.vertex)

  order <- c(in.v,out.v,inout.v)

  plot <- idealdigraph(wimp, ...)
    visIgraphLayout(plot,layout = "layout_as_star", center = center, order = order)
}
