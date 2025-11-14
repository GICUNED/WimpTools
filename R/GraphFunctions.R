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
#' @param vertex.vector Vector defining the value of each vertex in the digraph.
#'        Default is the normalized self from the wimp object.
#' @param ideal.vector Vector defining the ideal value of each vertex in the digraph.
#'        Default is the normalized ideal self from the wimp object.
#' @param width Digraph width.
#' @param height Digraph height.
#' @param color Color palette to be used. The options are "red/green" and "grey scale". Default is "red/green".
#' @param layout Layout with which the digraph will be displayed. The options
#'        are: "circle", "rtcircle", "tree", "graphopt", "mds" and "grid". Default is "graphopt".
#' @param show Logical vector defining which constructs to display. By default all constructs are shown.
#' @param hide.direct If TRUE, hide direct relationship between nodes of the graph. Default is FALSE.
#' @param areas Logical; if TRUE, draw coloured areas that group nodes by an attribute. Default is FALSE.
#' @param area.attr Character; the name of the vertex attribute to use for grouping nodes into areas.
#'        Default is "category". Can be any attribute column from the imported wimp vertices.
#' @param pad.side Numeric; padding in pixels used for the areas drawn around the nodes. Default is 50.
#' @param rounding Numeric; radius used to round the area corners (0 for sharp corners). Default is 10.
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
                    height="700px", color = "red/green", layout ="graphopt",
                    show = TRUE, hide.direct = FALSE,
                    areas = FALSE, area.attr = "category", area.color = NA, pad.side = 50, rounding = 10){

  # Helper: place nodes within a cell (internal)
  .place_in_cell <- function(ids, cx, cy, cell_size){
    n <- length(ids)
    if (n == 1L) {
      data.frame(id = ids, x = cx, y = cy)
    } else if (n == 2L) {
      r <- cell_size * 0.22
      ang <- c(pi/6, 5*pi/6)
      data.frame(id = ids,
                 x = cx + r * cos(ang),
                 y = cy + r * sin(ang))
    } else {
      r <- cell_size * 0.28
      ang <- seq(0, 2*pi, length.out = n+1L)[- (n+1L)]
      data.frame(id = ids,
                 x = cx + r * cos(ang),
                 y = cy + r * sin(ang))
    }
  }

  # Read wimp object
  if(inherits(wimp, "wimp")){
    stopifnot(!is.null(wimp$vertices), is.data.frame(wimp$vertices))
    lpoles <- wimp$vertices$lpole
    rpoles <- wimp$vertices$rpole
    poles <- paste(lpoles, "-", rpoles)
    wmatrix <- if(!is.null(wimp$global$wmatrix)) wimp$global$wmatrix else matrix(0, nrow(wimp$vertices), nrow(wimp$vertices))
    self <- wimp$vertices$self
    ideal <- wimp$vertices$ideal
    area_vec <- if(area.attr %in% names(wimp$vertices)) wimp$vertices[[area.attr]] else rep(NA_character_, nrow(wimp$vertices))
  } else if(inherits(wimp, "scn")){
    lpoles <- wimp$constructs[[1]]
    rpoles <- wimp$constructs[[2]]
    poles <- paste(lpoles, "-", rpoles)
    wmatrix <- wimp$weights
    self <- wimp[["self"]][[1]]
    ideal <- wimp[["self"]][[2]]
    area_vec <- rep(NA_character_, length(lpoles))
  }

  if(!exists("area_vec")) area_vec <- rep(NA_character_, length(lpoles))

  # Validate that area.attr exists; warn if not
  attr_exists <- if(!is.null(wimp$vertices) && is.data.frame(wimp$vertices)) area.attr %in% names(wimp$vertices) else FALSE
  if(!attr_exists && areas) {
    warning("Vertex attribute '", area.attr, "' not found. Areas layout will not group by this attribute.")
  }

  if(is.na(vertex.vector[1])) vertex.vector <-  self
  vertex.vector <- sapply(vertex.vector,.thr)

  if(is.na(ideal.vector[1])) ideal.vector <- ideal

  # Build node labels based on sign of self value
  vertex.name <- character(length(vertex.vector))
  for (n in seq_along(vertex.vector)) {
    x <- vertex.vector[n]
    vertex.name[n] <- if (x < 0) lpoles[n] else if (x > 0) rpoles[n] else poles[n]
  }

  # Compute node colors and groups
  congruency.vector <- vertex.vector/ideal.vector
  discrepant.color <- .color.selection(color)[1]
  congruent.color  <- .color.selection(color)[2]
  undefined.color  <- .color.selection(color)[3]
  dilemmatic.color <- .color.selection(color)[4]

  vertex.color <- sapply(congruency.vector, function(x) dplyr::case_when(
    x < 0 && x != -Inf ~ discrepant.color,
    x > 0 && x !=  Inf ~ congruent.color,
    x == 0            ~ undefined.color,
    is.infinite(x)    ~ dilemmatic.color,
    .default          =  dilemmatic.color))

  vertex.group <- sapply(congruency.vector, function(x) dplyr::case_when(
    x < 0 && x != -Inf ~ "Discrepant",
    x > 0 && x !=  Inf ~ "Congruent",
    x == 0            ~ "Undefined Self",
    is.infinite(x)    ~ "Dilemmatic",
    .default          =  "Dilemmatic"))

  vertex.vadjust <- sapply(abs(vertex.vector), function(x) -4*x^2-25*x-42)

  vertex <- data.frame(
    id = 1:length(vertex.vector),
    label = vertex.name,
    group = vertex.group,
    category = as.character(area_vec),
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

  # Ensure the chosen area attribute is present as a column for grouping
  vertex[[area.attr]] <- as.character(area_vec)

  # Reorient weight matrix by the sign of each node
  for (n in seq_along(vertex.vector)) {
    i <- vertex.vector[n]
    if(i != 0){
      s <- i/abs(i)
      wmatrix[,n] <- wmatrix[,n] * s
      wmatrix[n,] <- wmatrix[n,] * s
    }
  }

  if(hide.direct){
    logical.dilemmatic <- ideal == 0
    wmatrix[wmatrix > 0] <- 0
    wmatrix[logical.dilemmatic,] <- 0
    wmatrix[,logical.dilemmatic] <- 0
  }

  # Build edges table from weight matrix
  weight.vector <- c(); from.vector <- c(); to.vector <- c()
  for (i in 1:nrow(wmatrix)) for (j in 1:ncol(wmatrix)) {
    w <- wmatrix[i, j]
    if (w != 0) { weight.vector <- c(weight.vector,w); from.vector <- c(from.vector,i); to.vector <- c(to.vector,j) }
  }

  if(color!= "grey scale" ){
    edges.color  <- ifelse(weight.vector > 0 , "grey","#CD5C5C")
    edges.dashes <- FALSE
  } else {
    edges.color  <- "grey"
    edges.dashes <- weight.vector <= 0
  }

  edge.curved <- logical(length(weight.vector)); n <- 1
  for (N in 1:nrow(wmatrix)) for (M in 1:ncol(wmatrix)) {
    if(wmatrix[M,N] != 0 && wmatrix[N,M] != 0) edge.curved[n] <- TRUE
    if(wmatrix[N,M] != 0) n <- n + 1
  }


  edges <- data.frame(
    from = from.vector, to = to.vector,
    arrows = "to",
    dashes = edges.dashes,
    smooth = edge.curved,
    color = edges.color,
    title = round(weight.vector, 2)
  )

  # -------- LAYOUTS --------
  if (layout == "graphopt") layout <- "layout_with_graphopt"
  if (layout == "circle")   layout <- "layout_in_circle"
  if (layout == "tree")     layout <- "layout_as_tree"
  if (layout == "mds")      layout <- "layout_with_mds"
  if (layout == "grid")     layout <- "layout_on_grid"

  # Special layout: "areas" (grid placement by attribute)
  use_areas_layout <- identical(layout, "areas")
  if (use_areas_layout) {
    cat2 <- vertex[[area.attr]]
    cat2[is.na(cat2) | cat2 == ""] <- "Sin categoría"
    vertex[[area.attr]] <- cat2
    cats <- unique(cat2); k <- length(cats)
    # Internal constants for grid placement
    ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k / ncol)
    cell <- 400

    centers <- lapply(seq_len(k), function(i){
      r <- floor((i-1)/ncol); c <- (i-1) %% ncol
      c(cx = (c - (ncol-1)/2) * cell,
        cy = ((nrow-1)/2 - r) * cell)
    })
    centers <- do.call(rbind, centers); rownames(centers) <- cats

    pos_list <- lapply(cats, function(cat){
      ids <- vertex$id[vertex[[area.attr]] == cat & !vertex$hidden]
      if (length(ids) == 0L) return(NULL)
      cx <- centers[cat, "cx"]; cy <- centers[cat, "cy"]
      .place_in_cell(ids, cx, cy, cell)
    })
    pos_df <- do.call(rbind, pos_list)

    vertex$x <- NA_real_; vertex$y <- NA_real_
    if (!is.null(pos_df)) {
      vertex$x[match(pos_df$id, vertex$id)] <- pos_df$x
      vertex$y[match(pos_df$id, vertex$id)] <- pos_df$y
    }
    # Remaining nodes: jitter near the origin
    missing <- is.na(vertex$x)
    if (any(missing)) {
      set.seed(33)
      vertex$x[missing] <- rnorm(sum(missing), 0, cell*0.05)
      vertex$y[missing] <- rnorm(sum(missing), 0, cell*0.05)
    }
  }

  # Render network
  if (use_areas_layout) {
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 0, labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = TRUE, multiselect = TRUE) %>%
      visPhysics(enabled = FALSE)
  } else if (layout == "rtcircle") {
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = "layout_as_tree", circular = TRUE) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 0, labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
  } else {
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = layout, randomSeed = 33) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 0, labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
  }

  # Ensure the viewport auto-centers using stabilized/afterDrawing events
  if (!use_areas_layout || !areas) {
    g <- g %>% visEvents(
      stabilized = htmlwidgets::JS(
        "function(){ if(!this._autoFitDone){ try{ this.fit(); }catch(e){} this._autoFitDone=true; } }"
      ),
      afterDrawing = htmlwidgets::JS(
        "function(ctx){ if(!this._autoFitDone){ try{ this.fit(); }catch(e){} this._autoFitDone=true; } }"
      )
    )
  } else {
    # For areas layout with fixed positions, fit once after stabilization
    g <- g %>% visEvents(
      stabilized = htmlwidgets::JS(
        "function(){ if(!this._autoFitDone){ try{ this.fit({animation: false}); }catch(e){} this._autoFitDone=true; } }"
      )
    )
  }

  # Ensure 'category' is kept in widget nodes
  if ("category" %in% names(vertex) && !is.null(g$x$nodes)) {
    if (NROW(g$x$nodes) == nrow(vertex)) g$x$nodes$category <- vertex$category
  }

  # Area overlays (convex hulls with label)
  if(areas){
    uniq_cats <- unique(stats::na.omit(vertex$category[vertex$category != ""]))
    if(length(uniq_cats) > 0){
      if (missing(area.color) || is.null(area.color) || all(is.na(area.color))) {
        base_cols <- c("#6EA8FE","#72D6A0","#F7B267","#D985B9","#8FD3FE",
                       "#B5E48C","#FFD166","#CDB4DB","#A0C4FF","#FFAFCC")
      } else {
        base_cols <- area.color
      }
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
        "function(ctx){
  var net=this;
  var PAD=", pad.side, ";
  var ROUND=", rounding, ";
  var LINE_W=2;

  var CAT = ", jsonlite::toJSON(uniq_cats, auto_unbox=TRUE), ";
  var FILL= ", jsonlite::toJSON(as.character(fills), auto_unbox=TRUE), ";
  var STROK=", jsonlite::toJSON(as.character(strokes), auto_unbox=TRUE), ";

  function isCCW(pts){var s=0;for(var i=0;i<pts.length;i++){var a=pts[i],b=pts[(i+1)%pts.length];s+=a.x*b.y-a.y*b.x;}return s>0;}
  function lineIntersect(p1,d1,p2,d2){var det=d1.x*d2.y-d1.y*d2.x;if(Math.abs(det)<1e-9) return null;var t=((p2.x-p1.x)*d2.y-(p2.y-p1.y)*d2.x)/det;return {x:p1.x+t*d1.x,y:p1.y+t*d1.y};}
  function convexHull(points){
    if(points.length<=1) return points.slice();
    var pts=points.slice().sort(function(a,b){return a.x!==b.x?a.x-b.x:a.y-b.y;});
    function cross(o,a,b){return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);}
    var lower=[],upper=[];
    for(var i=0;i<pts.length;i++){while(lower.length>=2 && cross(lower[lower.length-2], lower[lower.length-1], pts[i])<=0) lower.pop(); lower.push(pts[i]);}
    for(var j=pts.length-1;j>=0;j--){while(upper.length>=2 && cross(upper[upper.length-2], upper[upper.length-1], pts[j])<=0) upper.pop(); upper.push(pts[j]);}
    upper.pop(); lower.pop(); return lower.concat(upper);
  }
  function offsetConvex(pts, d){
    if(pts.length<3) return null;
    var H=pts.slice(); if(!isCCW(H)) H.reverse();
    var N=H.length, out=new Array(N);
    for(var i=0;i<N;i++){
      var p0=H[(i-1+N)%N], p1=H[i], p2=H[(i+1)%N];
      var e0={x:p1.x-p0.x,y:p1.y-p0.y}, e1={x:p2.x-p1.x,y:p2.y-p1.y};
      var l0=Math.hypot(e0.x,e0.y)||1, l1=Math.hypot(e1.x,e1.y)||1; e0.x/=l0; e0.y/=l0; e1.x/=l1; e1.y/=l1;
      var n0={x:e0.y,y:-e0.x}, n1={x:e1.y,y:-e1.x};
      var pB={x:p1.x+n0.x*d,y:p1.y+n0.y*d}, pC={x:p1.x+n1.x*d,y:p1.y+n1.y*d};
      var q=lineIntersect(pB,e0,pC,e1);
      out[i]= q ? q : {x:p1.x+(n0.x+n1.x)*d, y:p1.y+(n0.y+n1.y)*d};
    }
    return out;
  }
  function draw(ctx, pts, fill, stroke){
    if(!pts || pts.length<3) return;
    ctx.save(); ctx.globalCompositeOperation='source-over';
    ctx.fillStyle=fill; ctx.strokeStyle=stroke; ctx.lineWidth=LINE_W;
    if(ROUND<=0){
      ctx.beginPath(); ctx.moveTo(pts[0].x,pts[0].y);
      for(var i=1;i<pts.length;i++) ctx.lineTo(pts[i].x,pts[i].y);
      ctx.closePath(); ctx.fill(); ctx.stroke();
    }else{
      ctx.lineJoin='round'; ctx.lineCap='round'; ctx.miterLimit=2;
      ctx.beginPath(); var n=pts.length;
      for(var i=0;i<n;i++){
        var p0=pts[(i-1+n)%n], p1=pts[i], p2=pts[(i+1)%n];
        var v1x=p1.x-p0.x,v1y=p1.y-p0.y,v2x=p2.x-p1.x,v2y=p2.y-p1.y;
        var l1=Math.hypot(v1x,v1y)||1,l2=Math.hypot(v2x,v2y)||1;
        var rr=Math.min(ROUND,0.45*l1,0.45*l2);
        v1x/=l1; v1y/=l1; v2x/=l2; v2y/=l2;
        var p1_in={x:p1.x-v1x*rr,y:p1.y-v1y*rr}, p1_out={x:p1.x+v2x*rr,y:p1.y+v2y*rr};
        if(i===0) ctx.moveTo(p1_in.x,p1_in.y); else ctx.lineTo(p1_in.x,p1_in.y);
        ctx.arcTo(p1.x,p1.y,p1_out.x,p1_out.y,rr);
      }
      ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    ctx.restore();
  }
  function drawLabelAbove(ctx, pts, text, stroke){
    var xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    var minX = Math.min.apply(null,xs), maxX = Math.max.apply(null,xs);
    var minY = Math.min.apply(null,ys);
    var cx = (minX + maxX) / 2, y = minY - 12;
    ctx.save(); ctx.font='bold 14px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='bottom';
    ctx.strokeStyle = stroke; ctx.lineWidth = 4; ctx.strokeText(text, cx, y);
    ctx.fillStyle = 'rgba(0,0,0,0.90)'; ctx.fillText(text, cx, y); ctx.restore();
  }

  var nodes=net.body.data.nodes.get(); var bycat={};
  nodes.forEach(function(n){
    if(n.hidden) return;
    if(n.category==null) return;
    var c=String(n.category);
    if(c===''||c==='NA') return;
    (bycat[c]||(bycat[c]=[])).push(n.id);
  });

  for(var i=0;i<CAT.length;i++){
    var cat = CAT[i], ids = bycat[cat];
    if(!ids || ids.length===0) continue;
    var pos = net.getPositions(ids);
    var pts = ids.map(function(id){ return {x:pos[id].x, y:pos[id].y}; });

    var outline;
    if(pts.length<=2){
      var xs=pts.map(p => p.x), ys=pts.map(p => p.y);
      var minX=Math.min.apply(null,xs)-PAD, maxX=Math.max.apply(null,xs)+PAD;
      var minY=Math.min.apply(null,ys)-PAD, maxY=Math.max.apply(null,ys)+PAD;
      outline=[{x:minX,y:minY},{x:maxX,y:minY},{x:maxX,y:maxY},{x:minX,y:maxY}];
    } else {
      var hull=convexHull(pts); outline=offsetConvex(hull, PAD);
    }

    draw(ctx, outline, FILL[i], STROK[i]);
    drawLabelAbove(ctx, outline, String(cat), STROK[i]);
  }
}"
      )

      g <- g %>% visEvents(afterDrawing = htmlwidgets::JS(js_hulls))
    }
  }

  g
}

# Ideal Digraph -----------------------------------------------------------

#' Ideal digraph -- idealdigraph()
#'
#' @description Plot the ideal self based on the constructs and their relations.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param inc If TRUE, hide direct relationships between nodes. Default is FALSE.
#' @param ... Additional arguments passed to \code{\link{digraph}}.
#' @param layout Layout for the digraph. Options: "circle", "rtcircle", "tree",
#'        "graphopt", "mds", "grid" or "areas". Default is "circle".
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A visNetwork graph
#'
#' @export
#'
#' @examples
#' idealdigraph(example.wimp)

idealdigraph <- function(wimp, inc=FALSE, layout = "circle", ...){

  # Extract ideal vector from new-format vertices
  ideal.vector <- wimp$vertices$ideal
  plot <- digraph(wimp = wimp, hide.direct = inc, vertex.vector = ideal.vector, layout = layout, ...)
  return(plot)
}

# Simulation Digraph ---------------------------------------------------------

#' Simulation digraph -- simdigraph()
#'
#' @description Plot the hypothetical self for a given scenario iteration.
#'
#' @param scn A scenario matrix ("scn" S3 object) from \code{\link{scenariomatrix}}.
#' @param niter Iteration index to display (0-based). Default is 0.
#' @param ... Additional arguments passed to \code{\link{digraph}}.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A visNetwork graph
#'
#' @export
#'
#' @examples
#' scn <- scenariomatrix(example.wimp, rep(1,5))
#' simdigraph(scn, niter = 2)

simdigraph <- function(scn, niter = 0, ...){

  # Select the vertex vector for the requested iteration
  vertex.vector <- scn[[1]][niter+1,]
  digraph(wimp = scn, vertex.vector = vertex.vector, ...)
}

# In - Out Digraph -------------------------------------------------------------

#' In-Out digraph -- inout_digraph()
#'
#' @description Show in, out, and in-out vertices for a given construct.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param iso.con Construct index (integer) to analyze.
#' @param ... Additional arguments passed to \code{\link{digraph}}.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A visNetwork graph
#'
#' @export
#'
#' @examples
#' inout_digraph(example.wimp, iso.con = 1)

inout_digraph <- function(wimp, iso.con, ...){

  # Isolate the target construct and rebuild a smaller wimp
  wimp <- .isolate.construct(wimp, iso.con)

  # Extract absolute weight matrix
  if(!is.null(wimp$global$wmatrix)){
    wmatrix <- abs(wimp$global$wmatrix)
  } else if(!is.null(wimp$vertices) && is.data.frame(wimp$vertices)){
    # Reconstruct from edges if needed
    n <- nrow(wimp$vertices)
    wmatrix <- matrix(0, n, n)
    if(!is.null(wimp$edges)){
      for(r in seq_len(nrow(wimp$edges))){
        i <- wimp$edges[r, "from"]
        j <- wimp$edges[r, "to"]
        wmatrix[i, j] <- wimp$edges[r, "weight"]
      }
    }
    wmatrix <- abs(wmatrix)
  }
  
  # Compute central node (max total degree)
  center <- which.max(colSums(wmatrix) + rowSums(wmatrix))

  out.vertex <- which(colSums(wmatrix)!=0)
  in.vertex <- which(rowSums(wmatrix)!=0)
  inout.v <- intersect(in.vertex,out.vertex)

  out.v <- setdiff(out.vertex,in.vertex)
  in.v <- setdiff(in.vertex,out.vertex)

  order <- c(in.v,out.v,inout.v)

  # Plot ideal digraph and apply star layout order
  plot <- idealdigraph(wimp, ...)
  visIgraphLayout(plot, layout = "layout_as_star", center = center, order = order)
}