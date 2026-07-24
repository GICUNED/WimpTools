#' Generate Tabbed HTML Widget for Psychlab
#'
#' @description
#' Creates a standalone HTML widget containing a two-tab interface for Psychlab integration.
#' The first tab displays the interactive digraph network, and the second tab displays 
#' the weight matrix heatmap.
#'
#' @param x A \code{wimp} object created by \code{\link{wimp}}.
#' @param ... Additional arguments passed to \code{digraph}.
#'
#' @return A `browsable` HTML object containing the tabbed interface.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
#'
#' @examples
#' \dontrun{
#' # Assuming 'w' is a valid wimp object
#' widget <- widget_digraph(w)
#' widget
#' }
widget_digraph <- function(x, ...) {
  if (!inherits(x, "wimp")) stop("Input must be a 'wimp' object.")
  
  # 1. Generate individual widgets
  g <- digraph(x, ...)
  h <- weight_heatmap(x)
  
  # 2. Define custom CSS for tabs
  css <- "
    body, html {
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
    }
    * {
      box-sizing: border-box;
    }
    .wt-tab-container {
      width: 100%;
      height: 100vh;
      display: flex;
      flex-direction: column;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif;
    }
    .wt-tab-header {
      overflow: hidden;
      border: 1px solid #eaeaea;
      background-color: #ffffff;
      display: flex;
      border-radius: 4px 4px 0 0;
    }
    .wt-tab-header button {
      background-color: inherit;
      border: none;
      outline: none;
      cursor: pointer;
      padding: 12px 24px;
      transition: 0.3s;
      font-size: 14px;
      font-weight: 600;
      color: #666666;
      flex-grow: 1;
    }
    .wt-tab-header button:hover {
      background-color: #f9f9f9;
    }
    .wt-tab-header button.active {
      background-color: #ffffff;
      color: #8cc63f;
      border-bottom: 2px solid #8cc63f;
    }
    .wt-tab-content {
      display: none;
      padding: 0;
      border: 1px solid #ccc;
      border-top: none;
      flex-grow: 1;
      height: 0;
      background-color: #ffffff;
      border-radius: 0 0 4px 4px;
    }
    .wt-tab-content > .html-widget {
      width: 100% !important;
      height: 100% !important;
      flex-grow: 1;
    }
  "
  
  # 3. Define vanilla JS for tab switching
  js <- "
    function openPsychlabTab(evt, tabName) {
      var i, tabcontent, tablinks;
      
      // Hide all tab content
      tabcontent = document.getElementsByClassName('wt-tab-content');
      for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = 'none';
      }
      
      // Remove 'active' class from all buttons
      tablinks = document.getElementsByClassName('wt-tab-btn');
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(' active', '');
      }
      
      // Show the current tab, and add an 'active' class to the button that opened the tab
      document.getElementById(tabName).style.display = 'flex';
      document.getElementById(tabName).style.flexDirection = 'column';
      evt.currentTarget.className += ' active';
      
      // Trigger resize event so htmlwidgets (plotly, visNetwork) render correctly in newly visible containers
      window.dispatchEvent(new Event('resize'));
    }
  "
  
  # 4. Construct UI using htmltools
  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    
    htmltools::tags$div(class = "wt-tab-container",
      htmltools::tags$div(class = "wt-tab-header",
        htmltools::tags$button(
          class = "wt-tab-btn active", 
          onclick = "openPsychlabTab(event, 'psychlab_tab_graph')", 
          "Self Digraph"
        ),
        htmltools::tags$button(
          class = "wt-tab-btn", 
          onclick = "openPsychlabTab(event, 'psychlab_tab_heatmap')", 
          "Weight Matrix"
        )
      ),
      
      htmltools::tags$div(id = "psychlab_tab_graph", class = "wt-tab-content", style = "display:flex; flex-direction:column;", 
        g
      ),
      htmltools::tags$div(id = "psychlab_tab_heatmap", class = "wt-tab-content", 
        h
      )
    )
  )
  
  # 5. Make it browsable so it renders in RStudio viewer or standalone HTML
  return(htmltools::browsable(ui))
}
