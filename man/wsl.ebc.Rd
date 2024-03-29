% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wsl_ebc.R
\name{wsl.ebc}
\alias{wsl.ebc}
\title{Correct environmental bias of species observations via environmental clustering}
\usage{
wsl.ebc(
  obs = NULL,
  ras = NULL,
  pportional = TRUE,
  plog = TRUE,
  nclust = 50,
  sp.specific = TRUE,
  sp.cor = 0.5,
  keep.bias = TRUE,
  filter = FALSE,
  path = NULL,
  ...
)
}
\arguments{
\item{obs}{Object of class matrix or data frame with three columns named "sp.id" (character),
"x" (numeric) and "y" (numeric). More than one observation per species must be referenced. If
only one species is in the data.frame, its sole observations will be resampled (equally or
proportionally) across the environmental clusters.}

\item{ras}{Object of class RasterBrick, RasterStack, list of RasterLayer of desired resolution
and extent. Used to generate the map of clusters needed to summarize the environmental space of
the study area.}

\item{pportional}{Logical. Should environmental stratification of observations be proportional
to the clusters' areas? If TRUE, EBCp applies.}

\item{plog}{Logical. Should EBCp apply with a logarithm? If TRUE, a stratification consensus
between EBCp and EBCe applies.}

\item{nclust}{Number of chosen clusters. Default is 50.}

\item{sp.specific}{Logical. Should EBC apply only for species whose environmental bias follows
the overall one (i.e. the number of original species observations per cluster is correlated
with that of the full dataset)?}

\item{sp.cor}{If sp.specific = TRUE, spearman's correlation tests are by default set to 0.5;
i.e. species with r < 0.5 are excluded from the function outputs.}

\item{keep.bias}{Default is TRUE. Strongly recommended to use when sp.cor = TRUE. Per species,
should the number of observations of the cluster in which the species was originally sampled
most often, be preserved? Said differently, for each species, should the cluster with the most
original species observations be as representative as the cluster with the most corrected
observations? If TRUE, after EBC applies, the number of observations per species in their
densest original cluster is set to that of the densest corrected cluster.}

\item{filter}{Logical. Should the observations be filtered according to 'ras' resolution?}

\item{path}{Path folder where the new species observation files should be saved.}

\item{...}{Additional arguments passed on to the ‘clara’ function (package 'cluster')}
}
\value{
The function returns in ‘path’, one text file of corrected observations (presences
or absences) per species. If the number of new EBC observations per species is too large,
sampling those randomly without replacement before model calibrations is advised.
}
\description{
Our wsl.ebc function corrects environmental bias from an observational dataset based on
environmental stratification / clustering. A map of n clusters is generated based on input
raster layers (in the article, climate predictors used in PPMs are employed as inputs).
Following random equal-stratified sampling design (EBCe; Hirzel & Guisan, 2002), for each
cluster, the number of observations per species (relative to all others) is artificially
rescaled to the total number of observations found in the cluster presenting the highest
observation density. Proportional-stratified sampling design (EBCp; Hirzel & Guisan, 2002)
adds a second step: for each cluster, the number of observations per species may additionally
be multiplied by the cluster's area (i.e. proportion of pixels in percentage relative to the
study area), or by its logarithm (default; consensus between EBCe and EBCp). Resulting output
indicates a new number of observations per cluster and species, that the function automatically
sub-samples with replacement over the original observational dataset. This function may be used
for presences and absences distinctively.
}
\examples{


}
\references{
Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., & Thuiller, W. (2021).
Novel methods to correct for observer and sampling bias in presence‐only species distribution
models. Global Ecology and Biogeography, 30(11), 2312-2325.
}
\author{
Yohann Chauvier
}
