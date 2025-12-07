# Helper functions for testing

# Function to load test data consistently
load_test_data <- function() {
  load(file.path(system.file(package = "WimpTools"), "data", "su_wimp.RData"))
  return(su_wimp)
}

# Function to check if an object has wimp structure
is_wimp_like <- function(obj) {
  if (!is.list(obj)) return(FALSE)
  if (!all(c("global", "vertices") %in% names(obj))) return(FALSE)
  if (!is.list(obj$global)) return(FALSE)
  if (!is.data.frame(obj$vertices)) return(FALSE)
  return(TRUE)
}

# Function to check plotly output
is_plotly_output <- function(obj) {
  inherits(obj, c("plotly", "htmlwidget"))
}