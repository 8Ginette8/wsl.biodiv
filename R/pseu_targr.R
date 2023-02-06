### =========================================================================
### sampling proportional to presence observation distribution
### =========================================================================
#' Sample pseudo-absences proportional to target group distribution 2.0
#'
#' This is essentially an improved version of prop.sampling. It uses the 'density'
#' function from the package 'spatstat' to create a density surface of the supplied
#' point pattern and samples pseudo-absences from this density-distribution. It
#' additionally requires a raster layer (env.layer) as input, typically containing
#' an environmental variable used for modelling to extract information about the
#' extent of the study area and the resolution of the environmental information.
#' Furthermore, a mask of class Spatial can be supplied to constrain the
#' area within which pseudoabsences should be sampled. If no mask is supplied,
#' pseudoabsences are sampled in the non-NA cells of the env.layer.
#'
#' @param points Matrix or data.frame with column names 'x' and 'y' assumed to
#' be on the same scale and metric distances (m, km,...)
#' @param nsamples Number of pseudoabsences to be generated
#' @param p.avoid Same type and conditions as "points" argument; avoid sampling
#' pseudo-absences in cells where input coordinates fall
#' @param env.layer A raster layer of desired resolution and extent (and area
#' of interest -> non_NA cells)
#' @param mask Optional. An object of class Spatial to mask for sampling area of interest
#' @param ... Additional arugements to the density.ppp function (package 'spatstat')
#' @return nsamples x 2 matrix with drawn psuedo-absences
#' @examples
#'
#' # Load files
#'
#' data(AlpineConvention_lonlat)
#' data(exrst)
#' data(xy_ppm)
#' mypoints=xy.ppm[,c("x","y")]
#' id.avoid=sample(1:nrow(mypoints),500,replace=TRUE)
#' xy.avoid=mypoints[id.avoid,]
#' xy.process=mypoints[-id.avoid,]
#'
#' # Run pseu.targ()
#'
#' dens1=pseu.targr(points=xy.process,nsamples=5000,p.avoid=xy.avoid,
#'                  env.layer=rst[[1]],mask=shp.lonlat,adjust=0.1)
#' dens2=pseu.targr(points=xy.process,nsamples=5000,p.avoid=xy.avoid,
#'                  env.layer=rst[[1]],mask=shp.lonlat,adjust=0.5)
#' dens3=pseu.targr(points=xy.process,nsamples=1000,p.avoid=NULL,
#'                  env.layer=rst[[1]],mask=NULL,adjust=0.1)
#' dens4=pseu.targr(points=xy.process,nsamples=1000,p.avoid=NULL,
#'                  env.layer=rst[[1]],mask=NULL,adjust=0.5)
#'
#' # Plot results
#' 
#' par(mfrow=c(2,2))
#' 
#' plot(rst[[1]])
#' points(dens1,pch=20,cex=0.3)
#' plot(rst[[1]])
#' points(dens2,pch=20,cex=0.3)
#' plot(rst[[1]])
#' points(dens3,pch=20,cex=0.3)
#' plot(rst[[1]])
#' points(dens4,pch=20,cex=0.3)
#' 
#' @author Philipp Brun, Yohann Chauvier
#' @export
pseu.targr=function(points,nsamples=1000,p.avoid=NULL,env.layer,mask=NULL,...){

  ### ------------------------
  ### check input data
  ### ------------------------

  if(ncol(points)!=2 || !all(colnames(points)%in%c("x","y"))){

    stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
  }

  if(!(class(env.layer)%in%c("RasterBrick","RasterStack","RasterLayer"))){
    stop("env.layer should be of class RasterBrick, RasterStack or RasterLayer!")
  }

  ### ------------------------
  ### Insert a mask argument to create NAs where we don't want the sampling to be done
  ### Necessary argument = a spatial object
  ### ------------------------

  if(!(class(mask)[1]%in%c("NULL","SpatialPolygons","SpatialPoints","SpatialLines",
                           "SpatialPolygonsDataFrame","SpatialPointsDataFrame","SpatialLinesDataFrame"))){
    stop("The mask should be of class SpatialPolygons, SpatialPoints or SpatialLines!")
  }

  if(!is.null(mask)){
    env.layer=mask(env.layer,mask)
  }

  ### ------------------------
  ### To avoid sampling P-Abs in cells where choosen coordinates fall
  ### Necessary argument = same type as "points"
  ### ------------------------

  if (!is.null(p.avoid)){
    if(ncol(p.avoid)!=2 || !all(colnames(p.avoid)%in%c("x","y"))){
      stop("Parameter 'p.avoid' should be a data.frame/matrix with two columns named x and y!")
    }
  }

  ### ------------------------
  ### Prepare point pattern object
  ### ------------------------

  # Define Point Pattern object to calculate
  xt=extent(env.layer)
  owi=owin(xrange=c(xt@xmin,xt@xmax),yrange=c(xt@ymin,xt@ymax))
  myppp=ppp(x=points[,"x"],y=points[,"y"],window = owi)

  ### ------------------------
  ### Generate 'im' object with density info
  ### ------------------------

  x=sort(unique(coordinates(env.layer)[,1]))
  y=sort(unique(coordinates(env.layer)[,2]))

  dens=density(myppp,xy=list(x=x,y=y),...)
  rdens=raster(env.layer)
  values(rdens)=dens$v[nrow(dens$v):1,]

  ### ------------------------
  ### Draw locations proportional to point density
  ### ------------------------

  vls=values(rdens)*values(raster::area(rdens))

  # Only sample points where env data coverage is complete
  na.locs=is.na(values(env.layer))
  vls[na.locs]=0

  # Do not sample where your choosen presences points are
  if (!is.null(p.avoid)){
    # Check where on the rasters they occur
    p.where=rasterize(p.avoid,env.layer, fun='count')
    # Probability of 0 sampling on these locations
    pres.locs=!is.na(values(p.where))
    vls[pres.locs]=0
  }

  # Replace NA's with zero probability
  if(any(is.na(vls)) || any(vls<0)){
    vls[which(is.na(vls) | vls<0)]=0
  }

  # Sample from density distributions
  vls=vls/sum(vls)
  pts=sample(1:length(vls),nsamples,prob=vls,replace=T)

  # Determine coordinates of samples
  pt.out=coordinates(rdens)[pts,]
  colnames(pt.out)=colnames(points)

  return(pt.out)

}
