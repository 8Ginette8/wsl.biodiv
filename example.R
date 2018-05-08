# ###########################################################################
# Get an overview over the data and figure out which traits make most sense
# to be used
#
# $Date: 2017-11-29
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# National Institute of Aquatic Resources (DTU-Aqua), 
# Technical University of Denmark, Charlottenlund
# 
# ###########################################################################

### =========================================================================
### Initialise system
### =========================================================================
# Start with house cleaning
rm(list = ls()); graphics.off()

### =========================================================================
### Get ID
### =========================================================================

setwd("R:/brunp/cluster_jobs/SDM_correlation/Src/Tools/wsl_biodiv/")

### =========================================================================
### Initialise system
### =========================================================================

# Load necessary Packages
library(mgcv)
library(dismo)

# source functions
scr=list.files("functions/",full.names=T)

for(i in 1:length(scr)){source(scr[i])}

### =========================================================================
### Prepare data
### =========================================================================

# Take anguilla data set from dismo package
data("Anguilla_train")

# Run wsl glm
vrs=c("SegSumT","USRainDays","USSlope")

env=Anguilla_train[,vrs]
form=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))

modi=wsl.glm(pa=Anguilla_train$Angaus,
             env_vars = env,
             taxon="Angaus",
             replicatetype="cv",
             reps=5,
             project="prototest",
             formula=form)


