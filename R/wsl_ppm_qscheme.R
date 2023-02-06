### =========================================================================
### define wsl.ppm.qscheme function
### =========================================================================
#' To prepare PPM general quadratic environment with input species points
#' 
#' Not to be called directly by the user
#'
#' @author Yohann Chauvier
#' @export
wsl.ppm.qscheme = function(data,area.win,quads)
{
  # Setting species data in ppp class
  sp.ppp = ppp(data[,1],data[,2],window=area.win,check=FALSE)

  # Apply a quadrature scheme needed for PPM
  dim=rev(dim(area.win))
  Q = quadscheme(data=sp.ppp, dummy=quads, method="grid", ntile=dim, npix=dim)
  return(Q)
}
