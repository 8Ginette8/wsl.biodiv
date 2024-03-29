% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wsl_samplePseuAbs.R
\name{wsl.samplePseuAbs}
\alias{wsl.samplePseuAbs}
\title{Sample pseudoabsences using various strategies}
\usage{
wsl.samplePseuAbs(
  n = 10000,
  env.stack,
  type = "geographic",
  add.strat = 0,
  pres = numeric(),
  taxon = character(),
  geodist_fact = 1,
  geores_fact = 20,
  template_dir = NA,
  geo_nrep = 7,
  target.group_dir = NA,
  env.strat_path = NA,
  rAll = TRUE,
  force_spat_thin = "no",
  limdist = NA,
  set_max_npres_to_nabs = TRUE
)
}
\arguments{
\item{n}{number of pseudoabsence points desired. Default is 10000.}

\item{env.stack}{RasterStack with environmental layers for sampling and
extraction}

\item{type}{Desired sampling strategy. Options are 'geographic', 'density',
'random', 'target.group', 'geo.strat', 'env.strat' and 'env.semi.strat'
(see details). Default is 'geographic'}

\item{add.strat}{Fraction between 0 and 1; should strategy be complemented
by a fraction of environmental strata.}

\item{pres}{SpatialPoints object. Location of presence points. Necessary for
'geographic' sampling strategy and the best adding point if the downstream
functions will be used.}

\item{taxon}{Character; name of taxon of interest to keep track of meta
information.}

\item{geodist_fact}{Factor to adjust spatial autocorrelation lengths: for
'geographic' pseudoabsence point patterns, values below 1 increase
autocorrelation length; values above 1 decrease it; for 'density' sampling
it is the other way around.}

\item{geores_fact}{Aggregation factor 'geographic' template. Larger values
save calculation time, but decrease resolution of sampling points.}

\item{template_dir}{Directory where template raster for 'geographic'
sampling should be saved in/loaded from. If NA, nothing will be saved; if
provided, template will be saved in/loaded from directory depending on
whether a file already exists.}

\item{geo_nrep}{number of replicates of geographic models. More will create
a smoother pattern but increase computation time.}

\item{target.group_dir}{Directory where xy files of target group taxa are
stored. Must be supplied if sampling stragegy is 'target.group', must contain
a column names 'x' and 'y' with coordinates in the same projection as other
spatial data}

\item{env.strat_path}{Directory where sample of environmental strata for
'env.strat' sampling should be saved in/loaded from. If NA, nothing will be
saved; if provided, environmental strata will be saved in/loaded from
directory depending on whether a file already exists.}

\item{rAll}{should all data be read into memory for computation of environmental
strata? this is faster but you may run into memory issues for large rasters.}

\item{force_spat_thin}{Should minimum distance be enforced between points? Options
are 'no', 'presences', 'absences', and 'both'. By default thinning is defined for
pseudoabsences from 'geographic', 'density', 'random', and 'geo.strat' methods
with minimum distance according to the resolution of the template raster.
'presences' takes the minimum distance criterion from the template raster over
to the 'presence' points; 'absences' takes the criterion over to 'env.strat',
'env.semi.strat', and 'target.group'; 'both' does it for both.}

\item{limdist}{The minimum distance accepted for spatial thinning. Units
should be km if the spatial data is projected, otherwise the units of the
coordinate reference system used. If no value is supplied, the maximum
distance between two cells of the template raster will be taken.}

\item{set_max_npres_to_nabs}{logical. Should the maximum number of presences be
equal to the number of pseudoabsences defined.}
}
\value{
an object of class 'wsl.pseudoabsences' that can be plotted and
passed on to wsl.flex
}
\description{
Flexlible function to sample pseudoabsences with various strategies and
store the results in a 'wsl.pseudoabsences' object  that can be
passed on in the 'wsl.biodiv' pipeline.
}
\details{
'geographic' samples pseudoabsences with a sampling probability.
inversely proportional to the geographic distance to presence observations.
'density' samples pseudoabsences proportional to the density of presence
observations. 'random' samples pseudoabsences randomly with a sampling
probability proportional the area of the cells. 'target.group' samples
pseudoabsences from the presences of the taxa of the target group, attempting
to correct for sampling bias. It depends on a directory with taxa defined by the
user as target group. 'geo.strat' samples pseudoabsences geographically
stratified either on a plane, or on a sphere depending on the projection of
the supplied env.stack. 'env.strat' samples pseudoabsences environmentally
stratified. Points are sampled form all realized combinations of environmental
conditions occuring in the environmental stack that have a minimal occurrence
frequency. By default environmental strata are calculated based on all
raster layers supplied. If a directory is supplied as 'env.strat_path', a
large sample of stratified points will be saved to speed up computations for
follow-up species. If environmental strata based on different predictors than
supplied are preferred 'env.strat_path'can be an .RData file from a
previous sampling of strata form different environmental predictors. 'env.semi.strat'
is similar to 'env.strat' but samples environmental strata proportional to the
logarithm of the area they cover.
}
\examples{
### =========================================================================
### Data preparation
### =========================================================================

# Predictors
bio=getData('worldclim',var='bio',lon=16, lat=48,res=.5)
bio=bio[[c(1,4,12)]]

# install.packages("rgbif")
library(rgbif)
# extract species
spn='Boletus aestivalis'
xt=as.vector(extent(bio))
baest <- occ_search(scientificName=spn,
                  hasCoordinate=TRUE,
                  decimalLongitude=paste0(xt[1],",",xt[3]),
                  decimalLatitude=paste0(xt[2],",",xt[4]))

pbaest=baest$data[,c('decimalLongitude','decimalLatitude')]
baest_spp=SpatialPoints(pbaest,proj4string = crs(bio))

# extract target group
targr <- occ_search(familyKey = 8789,
                   hasCoordinate=TRUE,
                   limit = 10000,
                   decimalLongitude=paste0(xt[1],",",xt[3]),
                   decimalLatitude=paste0(xt[2],",",xt[4]))

ptargr=as.matrix(targr$data[,c('decimalLongitude','decimalLatitude')])
colnames(ptargr)=c("x","y")

# create temporary directory for target.group info
tdir=paste0(tempdir(),"/trgr")
dir.create(tdir)
write.table(ptargr,file=paste0(tdir,"/targetxy.txt"),row.names = F)

# create temporary directory for template raster and env strata
strdir=paste0(tempdir(),"/str")
dir.create(strdir)

# Note that for these should not be temporary files for a real analysis.

### =========================================================================
### Sample pseudoabsences
### =========================================================================

# Geograhpic method with 20\% env strata
pseu.abs1=wsl.samplePseuAbs(type="geographic",
                          n=5000,
                          env.stack=bio,
                          pres=baest_spp,
                          add.strat=0.2,
                          template_dir=strdir,
                          env.strat_path=strdir,
                          geodist_fact=1,
                          geores_fact=3,
                          geo_nrep=7,
                          taxon=spn)

plot(pseu.abs1)

# Only geographic with longer autocorrelation length
pseu.abs2=wsl.samplePseuAbs(type="geographic",
                           n=5000,
                           env.stack=bio,
                           pres=baest_spp,
                           add.strat=0,
                           template_dir=strdir,
                           env.strat_path=strdir,
                           geodist_fact=.5,
                           geores_fact=3,
                           geo_nrep=7,
                           taxon=spn)

plot(pseu.abs2)

# Random and thin presences
pseu.abs3=wsl.samplePseuAbs(type="random",
                           n=5000,
                           env.stack=bio,
                           template_dir=strdir,
                           pres=baest_spp,
                           geores_fact=3,
                           add.strat=0,
                           taxon=spn,
                           force_spat_thin="presences")

plot(pseu.abs3)

# Geo.start
pseu.abs4=wsl.samplePseuAbs(type="geo.strat",
                           n=5000,
                           env.stack=bio,
                           template_dir=strdir,
                           pres=baest_spp,
                           geores_fact=3,
                           add.strat=0,
                           taxon=spn)

plot(pseu.abs4)

# Target group with 20\% env strat
pseu.abs5=wsl.samplePseuAbs(type="target.group",
                         n=5000,
                         env.stack=bio,
                         template_dir=strdir,
                         target.group_dir=tdir,
                         env.strat_path=strdir,
                         geores_fact=3,
                         pres=baest_spp,
                         add.strat=0.2,
                         taxon=spn,
                         force_spat_thin="both")

plot(pseu.abs5)

# Environmental semi-stratified
pseu_abs6=wsl.samplePseuAbs(n = 5000,
                           env.stack=bio,
                           type = "env.semi.strat",
                           add.strat = 0,
                           pres = baest_spp,
                           taxon = spn,
                           template_dir=strdir,
                           env.strat_path=strdir)

plot(pseu_abs6)

# Environmental semi-stratified with min dist
pseu_abs7=wsl.samplePseuAbs(n = 5000,
                           env.stack=bio,
                           type = "env.semi.strat",
                           add.strat = 0,
                           geores_fact=3,
                           pres = baest_spp,
                           taxon = spn,
                           template_dir=strdir,
                           env.strat_path=strdir)

plot(pseu_abs7)

# Density dependent
pseu_abs8=wsl.samplePseuAbs(n = 5000,
                           env.stack=bio,
                           type = "density",
                           add.strat = 0,
                           pres = baest_spp,
                           taxon = spn,
                           geores_fact=3,
                           template_dir=strdir,
                           env.strat_path=strdir)

plot(pseu_abs8)

### =========================================================================
### Fit SDMs
### =========================================================================

# Define model settings
vrs=names(bio)
form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
form.tree=as.formula(Presence ~ .)

modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=FALSE),
           multi("gbm",list(formula=form.tree,
                distribution = "bernoulli",
                interaction.depth = 1,
                shrinkage=.01,
                n.trees = 3500),"gbm-simple"),
           multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE),
           multi("randomForest",list(formula=form.tree,ntree=500,maxnodes=NULL),"waud1"))

# Fit models using wsl.flex function
library(gam)
library(gbm)
library(randomForest)
modi5=wsl.flex(x=pseu.abs5,
              replicatetype="block-cv",
              reps=3,
              strata=sample(1:3,nrow(pseu.abs5@env_vars),replace=TRUE),
              project="multitest",
              mod_args=modinp)

### =========================================================================
### Evaluate
### =========================================================================

# Evaluate the models with wsl.evaluate function
eval5<-wsl.evaluate(modi5,crit="maxTSS",prevalence_correction = TRUE)
summary(eval5)

### =========================================================================
### Predict
### =========================================================================

modi_pred=wsl.flex(x=pseu.abs5,
                  replicatetype="none",
                  reps=1,
                  project="test_pred",
                  mod_args=modinp[c(1:2)])

# Get thresholds
thr.5=get_thres(eval5)[1:2]

### Make some predictions
pred5=wsl.predict(modi_pred,predat=bio,thres = thr.5)

par(mfrow=c(1,2))
plot(pred5@predictions[[1]]$`glm-simple`,main="GLM")
plot(pred5@predictions[[1]]$`gbm-simple`,main="GBM")

}
\references{
Descombes, P., Chauvier, Y., Brun, P., Righetti, D., Wüest, R. O., Karger, D. N., ... &
Zimmermann, N. E. (2022). Strategies for sampling pseudo-absences for species distribution
models in complex mountainous terrain. bioRxiv, 2022-03.
}
\author{
Philipp Brun
}
