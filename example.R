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
library(MASS)
library(gbm)

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

### =========================================================================
### Check out wsl.glm
### =========================================================================

form=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))

# Try out wsl.glm funcion
modi1=wsl.glm(pa=Anguilla_train$Angaus,
             env_vars = env,
             taxon="Angaus",
             replicatetype="cv",
             reps=5,
             project="prototest",
             mod_tag="test-glm",
             formula=form,
             family="binomial",
             step=T)

# Try out custom summary function
summary(modi1)

# Access glm object of first replicate
summary(modi1@fits$`test-glm01`)

### =========================================================================
### Check out wsl.gam
### =========================================================================

form=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))

# Try out wsl.glm funcion
modi2=wsl.gam(pa=Anguilla_train$Angaus,
             env_vars = env,
             taxon="Angaus",
             replicatetype="splitsample",
             reps=3,
             project="prototest",
             mod_tag="test-gam",
             formula=form,
             family="binomial",
             step=F)

# Try out custom summary function
summary(modi2)

# Access glm object of first replicate
summary(modi2@fits$`test-gam01`)

### =========================================================================
### Check out wsl.gbm
### =========================================================================

# Try out wsl.glm funcion
modi3=wsl.gbm(pa=Anguilla_train$Angaus,
              env_vars = env,
              taxon="Angaus",
              replicatetype="none",
              reps=1,
              project="prototest",
              mod_tag="test-brt",
              formula=Presence ~ . ,
              distribution = "bernoulli",
              interaction.depth = 1,
              shrinkage=.01,
              n.trees = 3500)

# Try out custom summary function
summary(modi3)

# Access glm object of first replicate
summary(modi3@fits$`test-brt01`)

### =========================================================================
### Check out wsl.maxent
### =========================================================================

feat=c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")

# Try out wsl.glm funcion
modi4=wsl.maxent(pa=Anguilla_train$Angaus,
                 env_vars = env,
                 taxon="Angaus",
                 replicatetype="block-cv",
                 reps=3,
                 strata=sample(1:3,nrow(env),replace=T),
                 project="prototest",
                 mod_tag="test-mxe",
                 args=feat)

# Try out custom summary function
summary(modi4)

# Access glm object of first replicate
summary(modi4@fits$`test-mxe01`)

