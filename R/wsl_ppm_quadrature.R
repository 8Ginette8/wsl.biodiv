#' Poisson Point Process Models (PPPM) 'quadrature'
#'
#' Set up quadrature points (or "background points") necessary to run 'wsl.ppmO' and
#' 'wsl.ppmGlasso' fit functions. Those points apply a spatial scaling proportional to the study
#' area and estimate the maximised model log likelihood (see Renner 2013, Renner et al. 2015)
#'
#' @param mask Object of class 'RasterLayer' or 'RasterBrick' or 'RasterStack'. Mask is defined
#' according to non NAs values. Defines study area in pixels used to grid sample the quadrature points.
#' Sampling is done over raster centroids so the resolution of the mask defines the desired sampling.
#' @param area.win Object of class 'owin'. 'owin' output of 'wsl.ppm.window'
#' @param random Logical. Should quadrature points be generated randomly or according
#' to a regular mask ?
#' @param nQ To choose the number of random quadrature points in case 'random=TRUE'
#' @param lasso Logical. Form of the output. If TRUE, 'wsl.ppmGlasso' is used
#' @param env_vars Only when 'lasso=TRUE'. Object of class 'RasterBrick' or 'RasterStack'.
#' Set of predictors the user wants to use to fit the model
#' @return Object of class 'ppp' or 'wsl.quads'. Points associated to NAs env. values are removed
#' @author Yohann Chauvier
#' @examples
#' 
#' #### Load
#' 
#' data(AlpineConvention_lonlat)
#' data(exrst)
#' 
#' ### Define mask
#' 
#' maskR = mask(rst[[1]],shp.lonlat)
#' 
#' ### Run 'wsl.ppm.window' function
#' 
#' wind = wsl.ppm.window(mask = maskR,
#'                       val = 1,
#'                       owin = TRUE)
#' 
#' ### Define quadrature points for 'wsl.ppmGlasso'
#' 
#'    # Grid regular
#' quadG1 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = FALSE,
#'                         lasso = TRUE,
#'                         env_vars = rst)
#' 
#'    # Randomly
#' quadG2 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = TRUE,
#'                         nQ = 100000,
#'                         lasso = TRUE,
#'                         env_vars = rst)
#' 
#' ### Define quadrature points for 'wsl.ppmO'
#' 
#'    # Grid regular
#' quadO1 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = FALSE,
#'                         lasso = FALSE,
#'                         env_vars = NULL)
#' 
#'    # Randomly
#' quadO2 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = TRUE,
#'                         nQ = 100000,
#'                         lasso = FALSE,
#'                         env_vars = NULL)
#'
#' @export
### ==================================================================
### To create the quadrature points
### ==================================================================

wsl.quadrature = function(mask,area.win,random=FALSE,nQ=100000,lasso=TRUE,env_vars=NULL)
{
  # Extract non NAs values from mask
  mask.cells = which(!is.na(mask[]))

  # If random sampling
  if (random) {
    mask.sample = try(sample(mask.cells,nQ,replace=FALSE),silent=TRUE)
    if (class(mask.sample)%in%"try-error") {
      mask.sample = sample(mask.cells,nQ,replace=TRUE)
    }
    mask.cells = mask.sample
  }

  # Create our XY points defining our study area
  quad = xyFromCell(mask,mask.cells)

  # Creating our quadrature points
  quadS = ppp(quad[,1],quad[,2],window=area.win)

  # Other type of quadrature points in case of lasso
  if (lasso) {

    # Extract environmental values on our quadrature points
    qind = cellFromXY(env_vars,quad)
    qenv = extract(env_vars,qind)
    quads = as.data.frame(cbind(Presence=0,quad,qenv))
    quads = quads[complete.cases(quads),]

    # Keep infos
    names(quads)[4:ncol(quads)] = names(env_vars)

    # Return in a wsl.quads obejct
    quadS<-wsl.quads()
    quadS@coords = quads[,c("x","y")]
    quadS@Qenv = quads[,-c(2,3)]
  }
  return(quadS)
}

