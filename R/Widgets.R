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
  
  export_name <- deparse(substitute(x))
  
  # 1. Generate individual widgets
  g <- digraph(x, export_name = export_name, ...)
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
    :root {
      accent-color: #8cc63f !important;
    }
    input[type=\"checkbox\"], input[type=\"radio\"], input[type=\"range\"] {
      accent-color: #8cc63f !important;
    }
    select:focus { 
      border-color: #8cc63f !important; 
      outline: none; 
    }
    option:checked { 
      background-color: #8cc63f !important; 
      color: white !important; 
    }
    option:hover, option:focus, option:active { 
      background-color: #8cc63f !important; 
      color: white !important;
      box-shadow: 0 0 10px 100px #8cc63f inset !important;
    }
    ::selection { 
      background-color: #8cc63f !important; 
      color: white !important; 
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
      position: relative;
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
    
    function downloadHeatmap() {
      var plotEl = document.querySelector('#psychlab_tab_heatmap .js-plotly-plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: 'WIMP_EXPORT_NAME_Weight_Matrix'});
      }
    }
  "
  js <- gsub("WIMP_EXPORT_NAME", export_name, js)
  
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
        h,
        htmltools::HTML("
          <div style='position: absolute; bottom: 10px; right: 50px; z-index: 1000; background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='Export PNG' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"downloadHeatmap()\">
            <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line></svg>
          </div>
          
          <div style='position: absolute; bottom: 10px; right: 90px; z-index: 1000; background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='Info' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"document.getElementById('heatmap_info_modal').style.display='block';\">
            <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line></svg>
          </div>
          
          <div style='position: absolute; bottom: 10px; right: 10px; z-index: 1000; background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='Fullscreen' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"var el=this.closest('.wt-tab-container')||this.closest('.wt-tab-content'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}\">
            <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path></svg>
          </div>
          
          <div id='heatmap_info_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 80%; max-width: 400px; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
            <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
              <h3 style='margin:0; color:#444; font-size:16px;'>Weight Matrix</h3>
              <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('heatmap_info_modal').style.display='none';\">&times;</span>
            </div>
            <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.</p>
          </div>
        ")
      )
    )
  )
  
  # 5. Make it browsable so it renders in RStudio viewer or standalone HTML
  return(htmltools::browsable(ui))
}
