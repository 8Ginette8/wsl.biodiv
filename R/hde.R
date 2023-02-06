### =========================================================================
### define hde function
### =========================================================================
#' Hide annoying prints
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun
#' @export
hde=function(x){
  invisible(capture.output(x))
}
