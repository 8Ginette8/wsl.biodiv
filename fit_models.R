# ###########################################################################
# Definition of class pfit.2 and functions pfiz and summary
#
# $Date: 2014-06-22 V1
# $Date: 2018-03-14 V2

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
pfit.2<-setClass("pfit.2",slots=c(author="character", # Who created the object
                                  date="Date", # When was it created
                                  taxon="character", # Test data
                                  input="list", # Input data
                                  maxobs="numeric", # Maximum allowed observations
                                  replicatetype="character", # How should be replicated
                                  replicates="numeric", # Number of replicates
                                  tesdat="list", # Test data subset
                                  model.matrix="matrix", # Which models should be fitted
                                  fits="list")) # Model objects

### =========================================================================
### define pfit function
### =========================================================================

## Explanations (of non-obvious things):
# mm = model matrix: includes information which models should be fitted
#      and is mainly relevant for special designs (i.e. if not the same models
#      are required for the different data classes). It has the form
#          me    glm   gam   brt
#   simple  1     1     1     1
#  complex  1     1     0     1

#      If not provided invormation will be generated from models %*% classes


# make pfit
fiz<-function(pa,
              env_vars,
              taxon,
              xy=matirx(),
              models=c("Maxent","GAM","GLM","BRT"),
              complexity=c("simple","complex"),
              maxobs=NA,
              replicatetype="none",
              reps=1,
              mm=NA){
  
  #generate pfit object
  out<-pfit.2()
  
  #add Meta info
  out@author<-basename(Sys.getenv("USERPROFILE"))
  out@date<-Sys.Date()
  out@taxon<-taxon
  
  #add maximum observation information
  if(is.na(maxobs)==T){maxobs<-nrow(env_vars)}
  
  # check maxobs argument
  if(maxobs>nrow(env_vars)){stop(paste0(maxobs," Observations chosen - only",nrow(env_vars)," available!"))}
  
  # add maxobs argument
  out@maxobs<-maxobs
  
  #add input data
  out@input<-list(pa=pa,env=env_vars,xy=xy)

  #add number of replicates
  out@replicates<-reps
  
  #check for replicate type
  if(replicatetype%in%c("none","crossvalidate","splitsample")==F){stop("Non-existing replicatetype!")}
  if(replicatetype%in%c("crossvalidate","splitsample")==T && (is.na("reps") ==T || reps<2)){stop("Give reasonalbe number of replicates!")}
  if(replicatetype=="none" && reps>1){reps=out@replicates=1;warning("Replicate type is 'none' but multiple replicates were chosen - replicates are set to 1")}
  
  #add replicate type
  out@replicatetype<-replicatetype
  
  ### ------------------------
  ### generate Model Matrix
  ### ------------------------ 
  rn=c("simple","complex")
  cn=c("Maxent","GAM","GLM","BRT")
  
  if(length(mm)==1 && is.na(mm)==T){
    mm<-matrix(0,nrow=2,ncol=4)
    colnames(mm)<-cn
    rownames(mm)<-rn
    for (i in 1:nrow(mm)){
      for (j in 1:ncol(mm)){
        if(rownames(mm)[i]%in%complexity && colnames(mm)[j]%in%models){
          mm[i,j]<-1
        }
      }
    }
  } else if(all(dim(mm)==c(2,4)) && all(mm%in%c(0,1))){
    colnames(mm)<-rn
    rownames(mm)<-cn
    
  } else {stop("Model Matrix (mm) wrongly specified!")}
  
  print(mm)
  out@model.matrix<-mm
  
  ### --------------------------------------------
  ### Generate subsets according to replicate type
  ### --------------------------------------------

  obschoice<-list()
  testing<-list()
  if(replicatetype=="none"){
    obschoice[[1]]<-sample(1:nrow(env_vars),maxobs,replace = FALSE)
    
  } else if (replicatetype=="splitsample"){
    for (i in 1:reps){
      obschoice[[i]]<-sample(1:nrow(env_vars),size=round(maxobs*0.7),replace=F)
      testing[[i]]<-c(1:nrow(env_vars))[-which(1:nrow(env_vars)%in%obschoice[[i]])]
    }
    
  } else {
    # Currently block cross-validation
    
    if(length(xy)==1 && is.na(xy)==T){stop("Coordinates required for block cross-validation!!")}
    
    subs<-sample(1:nrow(env_vars),maxobs)
    qnts<-quantile(xy[,1],probs=0:5/5)
      
    for (i in 1:reps){
      obschoice[[i]]<-subs[which(xy[subs,1]<qnts[i] | xy[subs,1]>=qnts[i+1])]
      testing[[i]]<-subs[which(xy[subs,1]>=qnts[i] & xy[subs,1]<qnts[i+1])]
    }
  }

  ### ------------------------
  ### Fit models
  ### ------------------------ 
  
  #combine model functions
  fncs=list(Maxent=run.me,GLM=run.glm,GAM=run.gam,BRT=run.brt)
  
  fit_i<-list()
  # Loop over complexity
  for(i in 1:nrow(mm)){
    
    fit_j<-list()
    # Loop over model types
    for(j in 1:ncol(mm)){
      
      if(mm[i,j]==1){
        
        fit_k<-list()
        # Loop over repilcates
        for (k in 1:length(obschoice)){
          
          print(paste("Fitting", rn[i], cn[j] ,"- replicate:", k))
          
          #Prepare data
          y<-obschoice[[k]]
          input<-cbind(pa[y],env_vars[y,,drop=F])
          colnames(input)[1]<-"Presence"
          
          # Fit the model
          fit_k[[k]]<-fncs[cn[j]][[1]](dfra=input,type=rn[i])
        }
        fit_j[[j]]<-fit_k
        
      } else {
        fit_j[[j]]<-list()
      }
      
    }
    names(fit_j)<-cn
    fit_i[[i]]<-fit_j

  }
  names(fit_i)<-rn
  out@fits<-fit_i  

  ### ------------------------
  ### Generate testing data
  ### ------------------------
  
  if(replicatetype%in%c("crossvalidate","splitsample")){
    t.dat<-lapply(testing,function(y){
      weg<-cbind(pa[y],env_vars[y,])
      colnames(weg)[1]<-"Presence"
      return(weg)
    })
    out@tesdat<-t.dat
  }
  return(out)
}

### =========================================================================
### Summary for pfit
### =========================================================================

setMethod("summary",signature(object="pfit.2"),definition=function(object){
  mm2summ<-vector()
  for(i in 1:4){
    for(j in 1:2){
      if (object@model.matrix[j,i]==1){
        mm2summ<-append(mm2summ,paste0(rownames(object@model.matrix)[j],
                                      colnames(object@model.matrix)[i]))
      }
    }
  }
  
  summ<-c("pfit.2 object generated at "=as.character(object@date),
          "Author is "=object@author,
          "Taxon is "=object@taxon,
          "Environmental variables are "=paste(colnames(object@input$env),collapse=", "),
          "Replicate type is "=object@replicatetype,
          "Number of replicates is "=object@replicates,
          "Fitted models are "=paste(mm2summ,collapse=", "))
  # dimnames(summ)<-list(rep("",6),rep("",2))
  return(summ)
})