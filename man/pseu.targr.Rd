% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseu_targr.R
\name{pseu.targr}
\alias{pseu.targr}
\title{Sample pseudo-absences proportional to target group distribution 2.0}
\usage{
pseu.targr(
  points,
  nsamples = 1000,
  p.avoid = NULL,
  env.layer,
  mask = NULL,
  ...
)
}
\arguments{
\item{points}{Matrix or data.frame with column names 'x' and 'y' assumed to
be on the same scale and metric distances (m, km,...)}

\item{nsamples}{Number of pseudoabsences to be generated}

\item{p.avoid}{Same type and conditions as "points" argument; avoid sampling
pseudo-absences in cells where input coordinates fall}

\item{env.layer}{A raster layer of desired resolution and extent (and area
of interest -> non_NA cells)}

\item{mask}{Optional. An object of class Spatial to mask for sampling area of interest}

\item{...}{Additional arugements to the density.ppp function (package 'spatstat')}
}
\value{
nsamples x 2 matrix with drawn psuedo-absences
}
\description{
This is essentially an improved version of prop.sampling. It uses the 'density'
function from the package 'spatstat' to create a density surface of the supplied
point pattern and samples pseudo-absences from this density-distribution. It
additionally requires a raster layer (env.layer) as input, typically containing
an environmental variable used for modelling to extract information about the
extent of the study area and the resolution of the environmental information.
Furthermore, a mask of class Spatial can be supplied to constrain the
area within which pseudoabsences should be sampled. If no mask is supplied,
pseudoabsences are sampled in the non-NA cells of the env.layer.
}
\examples{

# Load files

data(AlpineConvention_lonlat)
data(exrst)
data(xy_ppm)
mypoints=xy.ppm[,c("x","y")]
id.avoid=sample(1:nrow(mypoints),500,replace=TRUE)
xy.avoid=mypoints[id.avoid,]
xy.process=mypoints[-id.avoid,]

# Run pseu.targ()

dens1=pseu.targr(points=xy.process,nsamples=5000,p.avoid=xy.avoid,
                 env.layer=rst[[1]],mask=shp.lonlat,adjust=0.1)
dens2=pseu.targr(points=xy.process,nsamples=5000,p.avoid=xy.avoid,
                 env.layer=rst[[1]],mask=shp.lonlat,adjust=0.5)
dens3=pseu.targr(points=xy.process,nsamples=1000,p.avoid=NULL,
                 env.layer=rst[[1]],mask=NULL,adjust=0.1)
dens4=pseu.targr(points=xy.process,nsamples=1000,p.avoid=NULL,
                 env.layer=rst[[1]],mask=NULL,adjust=0.5)

# Plot results

par(mfrow=c(2,2))

plot(rst[[1]])
points(dens1,pch=20,cex=0.3)
plot(rst[[1]])
points(dens2,pch=20,cex=0.3)
plot(rst[[1]])
points(dens3,pch=20,cex=0.3)
plot(rst[[1]])
points(dens4,pch=20,cex=0.3)

}
\author{
Philipp Brun, Yohann Chauvier
}
