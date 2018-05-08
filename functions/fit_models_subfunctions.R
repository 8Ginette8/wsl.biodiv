# ###########################################################################
# Second order functions for model fitting 
#
# $Date: 2018-05-08 V2

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 

# ###########################################################################


### =========================================================================
### retrieve meta information
### =========================================================================

meta.info=function(env = parent.frame()){
  
  e=as.list(env)
  
  
  m.i=list()
  
  m.i$author=Sys.info()[["user"]]
  m.i$date=Sys.time()
  m.i$replicatetype=e$replicatetype
  m.i$replicates=e$reps
  m.i$taxon=e$taxon
  print(class(e$env_vars))
  m.i$env_vars=colnames(e$env_vars)
  
  m.i$project=e$project
  
  return(m.i)
}

### =========================================================================
### check input data
### =========================================================================

checks=function(env=parent.frame()){
  
  
  #check for replicate type
  if(env$replicatetype%in%c("none","cv","block-cv","splitsample")==F){stop("Non-existing replicatetype!")}
  
  if(env$replicatetype=="block-cv" && is.na(env$strata)){stop("Stratum vector needed for block crossvalidation!")}
  
  if(env$replicatetype=="block-cv" && unique(env$strata)!=env$reps){stop("Stratum vector has wrong number of levels!")}
  
  if(env$replicatetype%in%c("cv","block-cv","splitsample")==T && (is.na("reps") ==T || env$reps<2)){stop("Give reasonalbe number of replicates!")}
  
  if(env$replicatetype=="none" && env$reps>1){env$reps=1;warning("Replicate type is 'none' but multiple replicates were chosen - replicates are set to 1")}
  
  
  #check for file path
  if(env$save && is.na(env$project)){stop("supply project in which data should be saved!")}
  
  if(env$save && !(env$project%in%list.dirs(env$path))){stop(paste("Project directory not existing in",env$path,
                                                                   "- please create manually"))}  
  
  
}

### =========================================================================
### Subset data according to replicate types
### =========================================================================

partition=function(replicatetype,reps,dat,strata){
  
  obschoice<-list()
  testing<-list()
  if(replicatetype=="none"){
    obschoice[[1]]<-dat
    
  } else if (replicatetype=="splitsample"){
    for (i in 1:reps){
      obschoice[[i]]<-dat[sample(1:nrow(dat),size=round(maxobs*0.7),replace=F),]
      testing[[i]]<-dat[c(1:nrow(dat))[-which(1:nrow(dat)%in%obschoice[[i]])],]
    }
    
  } else if (grepl("cv",replicatetype)){
    
    if(replicatetype=="cv"){
      unistr=sample(1:5,size=nrow(dat),replace=T)      
    } else {
      unistr=unique(strata)      
    }
    
    for (i in 1:reps){
      obschoice[[i]]<-dat[which(strata!=unistr[i]),]
      testing[[i]]<-dat[which(strata!=unistr[i]),]
    }
    
  }
  
  out=list(training=obschoice,testing=testing)
  
  return(out)
}

