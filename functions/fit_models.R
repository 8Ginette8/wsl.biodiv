# ###########################################################################
# Functions and classes related to model fitting
#
# $Date: 2018-05-08 V2

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
#

# ###########################################################################


### =========================================================================
### Set wsl.fit Class
### =========================================================================

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
                  save=FALSE,
                  project=NA,
                  path=NA,
                  step=FALSE,
                  mod_tag="",
                  ...){

  # check and prepare data and output
  lis=preps(call=match.call())
  
  #loop over replicates
  fits=list()
  for(i in 1:reps){
    modi=list()
    modi[[1]]=glm(...,data=lis$train[[i]])
    if(step){
      modi[[1]]=stepAIC(modi[[1]],direction="both",trace=FALSE)
    }
    names(modi)=ifelse(mod_tag=="","glm",mod_tag)
    fits[[i]]<-modi
    
  }
  
  # Name the fits
  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))
  
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
                  save=FALSE,
                  project=NA,
                  path=NA,
                  step=FALSE,
                  mod_tag="",
                  ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  #loop over replicates
  fits=list()
  for(i in 1:reps){
    
    modi=list()
    modi[[1]]=gam(...,data=lis$train[[i]])
    if(step){
      modi[[1]]=step(modi[[1]],direction="both",trace=FALSE)
    }
    names(modi)=ifelse(mod_tag=="","gam",mod_tag)
    fits[[i]]<-modi
  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))
  
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
                     save=FALSE,
                     project=NA,
                     path=NA,
                     mod_tag="",
                     ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  # loop over replicates
  fits=list()
  for(i in 1:reps){
    
    modi=list()
    
    # Create directory for temporary MaxEnt data
    hde(me.temp.dir<-paste("tmp_Maxent",print(as.numeric(Sys.time())*1000,digits=15),sep="_"))
    dir.create(me.temp.dir)
    
    d.in<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
    vec<-lis$train[[i]][,"Presence"]
    
    modi[[1]]=maxent(x=d.in,p=vec,...,path=me.temp.dir)
    
    names(modi)=ifelse(mod_tag=="","mxe",mod_tag)
    fits[[i]]<-modi
    
    # Remove Temporary folder for Maxent
    unlink(me.temp.dir,recursive=TRUE)
  }
  
  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
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
                  save=FALSE,
                  project=NA,
                  path=NA,
                  mod_tag="",
                  ...){
  
  # check and prepare data and output
  lis=preps(call=match.call())
  
  # loop over replicates
  fits=list()
  for(i in 1:reps){
    
    modi=list()
    modi[[1]]=gbm(...,data=lis$train[[i]])
    
    names(modi)=ifelse(mod_tag=="","gbm",mod_tag)
    fits[[i]]<-modi
  }
  
  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits=fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}

### =========================================================================
### define wsl.multi.fit function
### =========================================================================

wsl.multi.fit<-function(pa=numeric(), 
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  mod_args=list()){
  
  # Check supplied model types
  for(i in 1:length(mod_args)){
    if(!(mod_args[[i]]@mod%in%c("glm","gam","gbm","maxent","randomForest"))){
      warning(paste(mod_args[[i]]@mod,"not in focal model functions. You might run in to problems when evaluating/predicting..."))
    }
  }
  
  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){
    
    modi=list()
    # loop over models
    for(j in 1:length(mod_args)){
      
      if(mod_args[[j]]@mod=="maxent"){
        
        # Create directory for temporary MaxEnt data
        hde(mod_args[[j]]@args$me.temp.dir<-paste("tmp_Maxent",
                                                  print(as.numeric(Sys.time())*1000,digits=15),
                                                  sep="_"))
        dir.create(mod_args[[j]]@args$me.temp.dir)
        
        mod_args[[j]]@args$x<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
        mod_args[[j]]@args$p<-lis$train[[i]][,"Presence"]
        
        hde(modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args))
        
        #Remove Temporary folder for Maxent
        unlink(mod_args[[j]]@args$me.temp.dir,recursive=T)
        
      } else {
        
        mod_args[[j]]@args$data=lis$train[[i]]
        
        if(mod_args[[j]]@mod=="randomForest"){
          mod_args[[j]]@args$data$Presence=as.factor(mod_args[[j]]@args$data$Presence)
        }
        
        modi[[j]]=do.call(mod_args[[j]]@mod,mod_args[[j]]@args)
      }
      
      names(modi)[j]=ifelse(mod_args[[j]]@tag=="",paste0("model_",j),mod_args[[j]]@tag)

      if(mod_args[[j]]@step){
        modi[[j]]=stepAIC(modi[[j]],direction="both",trace=FALSE)
      }
    }
    
    fits[[i]]=modi

  }
  
  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))
  
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
  df=as.data.frame(object@meta[c("author","date","project")])  
  
  if(as.character(object@call)[1]=="wsl.multi.fit"){
    df$model_tags=paste(names(object@fits[[1]]),collapse=", ")
  } else {
    df$model_tags=object@meta["model_tag"]    
  }

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
