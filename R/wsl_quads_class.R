### =========================================================================
### Set wsl.quads Class
### =========================================================================
#' An S4 class to store quadrature  objects
#'
#' @slot coordinates informations on the generated quadratures
#' @slot predictors values extracted for each quadrature
#' @author Yohann Chauvier, Philipp Brun
#' @export
wsl.quads<-setClass("wsl.quads",slots=c(coords="list", # Coordinates information 
                                    Qenv="list")) # Q env infromations
