# ###########################################################################
# Functions and classes related to model fitting
#
# $Date: 2018-05-08 V2

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
# Description: Here commands relating to the class pfit.2 are defined
# This is an extension of the previous class pfit that allows a more flexible
# Choice of models and model classes
#

# ###########################################################################


### =========================================================================
### Set Class
### =========================================================================

# generate pfit class
wsl.fit<-setClass("wsl.fit",slots=c(meta="list", # Meta information
                                    tesdat="list", # Test data subset
                                    fits="list")) # Model objects

### =========================================================================
### define pfit function
### =========================================================================

wsl.glm<-function(pa=numeric(), 
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=F,
                  project=NA,
                  path=NA,
                  step=F,
                  ...){
  
  # check input
  checks()
  
  # generate wsl.fit object
  out<-wsl.fit()
  
  # supply meta info
  out@meta<-meta.info()
  
  print("yes..")
  # partition observations according to replicate type
  dat=cbind(data.frame(Presence=pa),env_var)
  parts=partition(replicatetype,reps,dat,strata)
  
  # supply testing data
  out@tesdat=parts$testing
  
  #fit glm
  fits=list()
  for(i in 1:reps){
    
    fits[[i]]=glm(...,data=parts$training[[i]])
    if(step){
      fits[[i]]=stepAIC(fits[[i]],direction="both",trace=F)
    }
  }
  
  # supply fitted objects
  out@fits=fits
  
  
  # Save
  #...
  
  return(out)
  
}
