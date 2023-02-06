#' Poisson Point Process Models (PPPM) 'quadrature' for presence-absence evaluation metrics
#'
#' Sample independent test quadrature points over the study area used for evaluating presence-only
#' models with presence-absence metrics (i.e. AUC, TSS, Kappa, PPV...). Presence-only models can
#' be evaluated with other metrics than Boyce index by using random sampling points as absences.
#' But these absences validation data need to be generated independenly over the study area i.e.
#' by avoiding resampling quadrature points already used in the models.
#'
#' @param mask Object of class 'RasterLayer' or 'RasterBrick' or 'RasterStack' of. Mask is defined
#' according to non NAs values. Defines study area in pixels used to grid sample the quadrature points.
#' Sampling is done over raster centroids so the resolution of the mask defines the desired sampling.
#' Same as in 'wsl.quadrature'
#' @param quadrature Object of class 'wsl.quads' or 'ppp'. May be one object or a list of object.
#' Same as in wsl.ppm.fit
#' @param env_vars An object of class 'RasterStack' or 'RasterBrick'. Same spatial layers used in
#' the models. Use to extract env values from test quadrature points
#' @param nQ Mumber of random quadrature points to sample
#' @param replace Logical. If TRUE, quadrature sampling is done with replacements
#' @return Object of class 'wsl.quads' or 'list'. Points associated to NAs env. values are removed
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
#' ### Define quadrature points for 'wsl.ppmGlasso'
#' 
#' quadG2 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = TRUE,
#'                         nQ = 10000,
#'                         lasso = TRUE,
#'                         env_vars = rst)
#' 
#' ### Define quadrature points for 'wsl.ppmO'
#' 
#' quadO2 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = TRUE,
#'                         nQ = 10000,
#'                         lasso = FALSE,
#'                         env_vars = NULL)
#' 
#' ### Define quadrature test points for 'wsl.ppmGlasso'
#' 
#' test.quadG2=wsl.quads.test(mask=maskR,
#'                            quadrature=quadG2,
#'                            nQ=1000,
#'                            replace=TRUE)
#' 
#' 
#' ### Define quadrature test points for 'wsl.ppmO'
#' 
#' test.quadG2=wsl.quads.test(mask=maskR,
#'                            quadrature=quadO2,
#'                            nQ=1000,
#'                            replace=TRUE)
#' 
#' 
#' @export
### ==================================================================
### To create the quadrature points
### ==================================================================

wsl.test.quads = function(mask,quadrature,env_vars,nQ=10000,replace=TRUE)
{
  QQ=quadrature
  if (class(QQ)%in%"wsl.quads"|class(QQ)%in%"ppp") {
    quadrature=list(quadrature)
  }

  quads.pa=list()
  for (o in 1:length(quadrature))
  {
    if (class(quadrature[[o]])%in%"wsl.quads") {
      quads=quadrature[[o]]@coords
    } else {
      quads=coords(quadrature[[o]])
    }
    
    MASK=mask

    # Transform in NAs from the ref grid locations where we find quadrature points
    qpa1=wsl.quads()
    cell.quads=cellFromXY(MASK,quads)
    MASK[][cell.quads]=NA

    # Extract randomly from this grid non NAs points (probability of 1)
    MASK[][!is.na(MASK[])]=1
    MASK[][is.na(MASK[])]=0
    quad.id=sample(1:length(MASK[]),size=nQ,prob=MASK[],replace=replace)

    # Store information in wsl.quads object
    qpa1.xy=xyFromCell(MASK,quad.id)
    qpa1.na=na.omit(data.frame(Presence=0,qpa1.xy,extract(env_vars,qpa1.xy)))
    qpa1@coords=list(qpa1.na[,c(2,3)])
    qpa1@Qenv=list(qpa1.na[,-c(2,3)])
    quads.pa[[o]]=qpa1
  }
  return(quads.pa)
}

