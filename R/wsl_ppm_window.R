#' Poisson Point Process Models (PPPM) 'window'
#'
#' Set up a study window of class 'owin' or 'im' necessary to run 'wsl.ppmO' fit function.
#'
#' @param mask Object of class 'RasterLayer' or 'RasterBrick' or 'RasterStack'. Mask is defined
#' according to non NAs values. Defined study area in pixels. Same as in 'wsl.quadrature'
#' @param val Values that are assigned to the 'owin' or 'im' mask
#' @param owin TRUE or FALSE
#' @return Object of class 'owin' or 'im'
#' @author Yohann Chauvier
#' @examples
#' 
#' ### Load
#' 
#' data(AlpineConvention_lonlat)
#' data(exrst)
#' 
#' ### Define mask
#' 
#' maskR = mask(rst[[1]],shp.lonlat)
#' 
#' ### Run function
#' 
#' wind = wsl.ppm.window(mask = maskR,
#'                       val = 1,
#'                       owin = TRUE)
#' 
#' @export
### ==================================================================
### To create object of class "owim" or "im"
### ==================================================================

wsl.ppm.window = function(mask,val=1,owin)
{ 
  # Extract values from raster and set "1" for non-NAs
  binValues = mask[]
  binValues[!is.na(binValues)] = val

  # Extract coordinates from raster extent
  yrow = unique(coordinates(mask)[,"y"])
  xcol = unique(coordinates(mask)[,"x"])

  # Create a new matrix from vector & reverse it to the south for compatibility purpose with im()
  inputM = matrix(binValues,dim(mask)[1],dim(mask)[2],byrow=T)
  flip.IT = apply(inputM,2,rev)

  # Create a new mask of class "im" + the quadrature points defined by our study area
  area.win = im(flip.IT, xcol=xcol, yrow=yrow)

  # Return object of class "owin" or "im"
  if (owin) {
    return(as.owin(area.win))
  } else {
    return(area.win)
  }
}

