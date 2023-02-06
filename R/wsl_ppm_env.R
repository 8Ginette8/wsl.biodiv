#' Poisson Point Process Models (PPPM) 'env'
#'
#' Convert envrionmental data for 'wsl.ppmO' fitting.
#'
#' @param env_vars Object of class 'RasterBrick' or 'RasterStack'. Same extent and resolution as 'mask'
#' @param mask Object of class 'RasterLayer' or 'RasterBrick' or 'RasterStack'. Mask is defined
#' according to non NAs values and the resolution of the study
#' @return List of objects of class 'im'
#' @author Yohann Chauvier
#' @examples
#' 
#' ### Load
#' 
#' data(exrst)
#' data(xy_ppm)
#' mypoints = xy.ppm[,c("x","y")]
#' 
#' ### Define mask
#' 
#' maskR = mask(rst[[1]],shp.lonlat)
#' 
#' ### Define your environments
#' 
#'    # For 'wsl.ppmGlasso' (observations focus)
#' envG = raster::extract(rst,mypoints)
#' 
#'    # For 'wsl.ppmO' (study area focus)
#' envO = wsl.ppm.env(rst,maskR)
#' 
#' @export
### ==================================================================
### To convert environmental data for ppm fitting
### ==================================================================

wsl.ppm.env = function(env_vars,mask)
{
  # Change mask values in "1" & apply it to the predictors
  names_data.env = names(env_vars)
  values(mask)[!is.na(values(mask))]=1
  mask.env = env_vars * mask
  names(mask.env) = names_data.env

  # Transform predictors in columns & remove "pixels" where we mask is NA
  data.env = mask.env[]
  ind_NA = is.na(mask[])
  if(any(ind_NA)) {data.env = data.env[!ind_NA,]}

  # Create for each predictor a spatial object 'im'
  var.list = lapply(1:ncol(data.env),
    function(x) wsl.ppm.window(mask=mask,val=data.env[,x],owin=FALSE))
  names(var.list) = colnames(data.env)
  return(var.list)
}