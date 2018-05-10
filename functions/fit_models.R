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
                                    fits="list", # Model objects
                                    call="call")) # conserve function call

### =========================================================================
### define wsl.glm function
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
                  mod_tag="",
                  ...){

  # check and prepare data and output
  lis=preps(call=match.call())
  
  #fit glm
  fits=list()
  for(i in 1:reps){
    
    fits[[i]]=glm(...,data=lis$train[[i]])
    if(step){
      fits[[i]]=stepAIC(fits[[i]],direction="both",trace=F)
    }
  }
  
  # Name the fits
  if(mod_tag==""){
    mod_tag="glm"
  }
  names(fits)=paste0(mod_tag,sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}

### =========================================================================
### define wsl.gam function
### =========================================================================

wsl.gam<-function(pa=numeric(), 
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=F,
                  project=NA,
                  path=NA,
                  step=F,
                  mod_tag="",
                  ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  #fit glm
  fits=list()
  for(i in 1:reps){
    
    fits[[i]]=gam(...,data=lis$train[[i]])
    if(step){
      fits[[i]]=stepAIC(fits[[i]],direction="both",trace=F)
    }
  }
  
  # Name the fits
  if(mod_tag==""){
    mod_tag="gam"
  }
  names(fits)=paste0(mod_tag,sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}

### =========================================================================
### define wsl.maxent function
### =========================================================================

wsl.maxent<-function(pa=numeric(), 
                     env_vars=data.frame(),
                     taxon=character(),
                     replicatetype=character(),
                     reps,
                     strata=NA,
                     save=F,
                     project=NA,
                     path=NA,
                     mod_tag="",
                     ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  # Create directory for temporary MaxEnt data
  id <- as.numeric(commandArgs(trailingOnly=TRUE))
  me.temp.dir<-paste("tmp_Maxent",gsub("\\:","\\.",Sys.time()),sep="_")
  me.temp.dir<-paste(gsub(" ","_",me.temp.dir),id,sep="_")
  dir.create(me.temp.dir)
  
  #fit maxent
  fits=list()
  for(i in 1:reps){
    
    d.in<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
    vec<-lis$train[[i]][,"Presence"]
    
    fits[[i]]=maxent(x=d.in,p=vec,...,path=me.temp.dir)
  }
  
  # Name the fits
  if(mod_tag==""){
    mod_tag="mxe"
  }
  names(fits)=paste0(mod_tag,sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
  #Remove Temporary folder for Maxent
  unlink(me.temp.dir,recursive=T)
  
  # Save
  #...
  
  return(lis$wslfi)
  
}


### =========================================================================
### define wsl.gbm function
### =========================================================================

wsl.gbm<-function(pa=numeric(), 
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=F,
                  project=NA,
                  path=NA,
                  mod_tag="",
                  ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  #fit maxent
  fits=list()
  for(i in 1:reps){

    fits[[i]]=gbm(...,data=lis$train[[i]])
  }
  
  # Name the fits
  if(mod_tag==""){
    mod_tag="gbm"
  }
  names(fits)=paste0(mod_tag,sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}

### =========================================================================
### define summary function for wsl.fit objects
### =========================================================================

sm=setMethod("summary",signature(object="wsl.fit"),definition=function(object){
  
  cat("Call: \n")
  print(object@call)
  
  cat("\nMeta information: \n")
  
  df=as.data.frame(object@meta[c("author","date","project","model_tag")])
  rownames(df)=""
  print(df)
  
  cat("\nVariables used: \n")
  
  df=as.data.frame(object@meta[c("taxon","env_vars")])
  rownames(df)=""
  print(df)
  
  cat("\nResampling: \n")
  
  df=as.data.frame(object@meta[c("replicatetype","replicates")])
  rownames(df)=""
  print(df)
  
  cat("\nOther: \n")
  
  if(as.character(object@call)[1]=="wsl.glm"){
    df=as.data.frame(object@meta[c("step")])
    rownames(df)=""
    print(df)
  }
  

})
