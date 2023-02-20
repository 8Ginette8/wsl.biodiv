### =========================================================================
### Test spatial auto-correlation in model ouptuts
### =========================================================================

#' Test spatial auto-correlation in model ouptuts
#'
#' A handful wrapper of moran.I() from the R package 'ape'. Based on the
#' data coordinates & model residuals, the function run a moran test to check if
#' there is any spatial auto-correlation (p<0.05) in the model.
#'
#' @param xy Object of class data.frame. Coordinates of the data points. First
#' and secon columns must be x and y respectively.
#' @param resid Object of class vector. Model residuals.
#' @return moran.I() test outputs
#' @author Yohann Chauvier
#' @examples
#'
#' # Load
#' data(AlpineConvention_lonlat)
#' data(exrst)
#' rst = rst[[1:6]]
#' data(xy_ppm)
#' mypoints = xy.ppm[,c("x","y")]
#' 
#' # Define mask
#' maskR = mask(rst[[1]],shp.lonlat)
#' 
#' # Run 'wsl.ppm.window' function
#' wind = wsl.ppm.window(mask = maskR,
#'                     val = 1,
#'                     owin = TRUE)
#' 
#' # nDefine quadrature points for 'wsl.ppmGlasso'
#' quadG1 = wsl.quadrature(mask = maskR,
#'                       area.win = wind,
#'                       random = FALSE,
#'                       lasso = TRUE,
#'                       env_vars = rst)
#' 
#' # Define your environments
#' envG = raster::extract(rst,mypoints)
#' 
#'   # Simple PPP (poly = FALSE & lasso=FALSE) + block-cross validation
#' ppm.simple = wsl.ppmGlasso(pres = mypoints,
#'                      quadPoints = quadG1,
#'                      asurface = raster::area(shp.lonlat)/1000,
#'                      env_vars = envG,
#'                      taxon = "species_eg2",
#'                      replicatetype = "none",
#'                      reps = 1,
#'                      lasso = FALSE)

#' # Test for auto-correlation
#' moran.I(xy = rbind(mypoints,quadG1@coords),
#'         resid = ppm.simple@fits[[1]][[1]]$residuals)
#' 
#' @export
moranI.test=function(xy,resid) {

  # Compute spatial distances
  inv.dist = 1/as.matrix(dist(xy))
  diag(inv.dist) = 0
  inv.dist[inv.dist%in%Inf|inv.dist%in%-Inf] = 0

  # Test auto-correlation
  moran.all = Moran.I(resid,inv.dist)

  # Return
  return(moran.all)
}
