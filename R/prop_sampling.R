### =========================================================================
### sampling proportional to presence observation distribution
### =========================================================================
#' Sample pseudo-absences proportional to target group distribution
#'
#' Uses the 'density' function from the package spatstat to create a density
#' surface of the supplied point pattern and samples pseudo-absences from this
#' density-distribution.
#'
#' @param points matrix or data.frame with column names 'x' and 'y' assumed to
#' be on the same scale and metric distances (m, km,...)
#' @param nsamples number of pseudoabsences to be generated
#' @param res resolution of the denstiy grid from which pseudo-absences are drawn.
#' default (1) corresponds to 1000 cells on the x-axis
#' @param ... arguments passed on to the 'density' function from package spatstat.
#' Note in particluar the argument 'adjust' which controls the kernel size in the
#' density interpolation.
#' @return nsamples x 2 matrix with drawn psuedo-absences
#' @author Philipp Brun
#' @examples
#'
#' ### Load points
#' 
#' data(mypts)
#' 
#' ### Test proportional sampling
#' 
#' pseu.1=prop.sampling(points=my.pts,nsamples=10000,adjust=1,res=1)
#' pseu.2=prop.sampling(points=my.pts,nsamples=10000,adjust=.5,res=1)
#' pseu.3=prop.sampling(points=my.pts,nsamples=10000,adjust=.2,res=1)
#' pseu.4=prop.sampling(points=my.pts,nsamples=10000,adjust=.1,res=1)
#'
#' par(mfrow=c(2,2))
#'
#' plot(pseu.1,pch=16,col="red",cex=.5)
#' points(my.pts)
#' plot(pseu.2,pch=16,col="red",cex=.5)
#' points(my.pts)
#' plot(pseu.3,pch=16,col="red",cex=.5)
#' points(my.pts)
#' plot(pseu.4,pch=16,col="red",cex=.5)
#' points(my.pts)
#'
#' @export
prop.sampling=function(points,nsamples=1000,res=1,...){

  ### ------------------------
  ### check input data
  ### ------------------------

  if(ncol(points)!=2 || !all(colnames(points)%in%c("x","y"))){

    stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
  }

  ### ------------------------
  ### Prepare point pattern object
  ### ------------------------

  # get x range
  xrng=range(points$x)

  # Determine order of magnitude of 1% of x range
  buff=0.01*abs(diff(xrng))
  oom=-10:10
  myoom=oom[which.min(abs(log10(buff)-oom))]

  # Extend x and y range by that magnitude to define observational window
  xrng[1]=floor(xrng[1]/10^myoom)*10^myoom
  xrng[2]=ceiling(xrng[2]/10^myoom)*10^myoom

  yrng=range(points$y)
  yrng[1]=floor(yrng[1]/10^myoom)*10^myoom
  yrng[2]=ceiling(yrng[2]/10^myoom)*10^myoom

  # Define Point Pattern object to calculate
  owi=owin(xrange=xrng,yrange=yrng)
  myppp=ppp(x=points[,"x"],y=points[,"y"],window = owi)

  ### ------------------------
  ### Generate 'im' object with density info
  ### ------------------------

  lo=res*c(1000,round(1000*(yrng[2]-yrng[1])/(xrng[2]-xrng[1])))

  x=seq(xrng[1],xrng[2],length.out = lo[1])
  y=seq(yrng[1],yrng[2],length.out = lo[2])

  dens=density(myppp,xy=list(x=x,y=y),...)

  ### ------------------------
  ### Draw locations proportional to point density
  ### ------------------------

  vls=as.vector(dens$v)

  # Replace NA's with zero probability
  if(any(is.na(vls)) || any(vls<0)){
    vls[which(is.na(vls) | vls<0)]=0
  }

  # Sample from density distributions
  pts=sample(1:length(vls),nsamples,prob=vls)

  # Determine coordinates of samples
  indi=arrayInd(pts,.dim=rev(lo))
  pt.out=cbind(x[indi[,2]],y[indi[,1]])
  colnames(pt.out)=colnames(points)

  return(pt.out)

}
