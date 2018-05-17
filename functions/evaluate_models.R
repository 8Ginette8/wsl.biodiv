# ###########################################################################
# Evaluation function and class definition
#
# $Date: 2018-05-08 V2

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
#
#############################################################################
### =========================================================================
### Set Class
### =========================================================================

wsl.evaluation<-setClass("wsl.evaluation",slots=c(meta="list", # Meta information
                                                  thres="numeric", # supply external threshold
                                                  performance="list")) # conserve function call

### =========================================================================
### define wsl.evaluate function
### =========================================================================
# 
wsl.evaluate<-function(x,tester=data.frame(),thres=numeric(),crit="pp=op"){
  
  ### ------------------------
  ### check tresholds
  ### ------------------------
  
  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
      }
    if(crit!="external"){
      warning("Assuming you want external tresholds to be used - setting crit='external'!")
      crit="external"
    }
  }  
  
  if(!(crit%in%c("pp=op","max","external"))){
    stop("Invalid threshold criterion chosen!")
  }
  
  
  ### ------------------------
  ### Check testing data and prepare for evaluation
  ### ------------------------ 
  
  if(x@meta$replicatetype=="none" && nrow(tester)==0){
    stop("External testing data must be supplied for replicatetype 'none'")
  } else if(x@meta$replicatetype%in%c("cv","block-cv","splitsample")) {
    outerloop<-length(x@tesdat)
    testa<-lapply(x@tesdat,function(x){
      y<-x[,-which(colnames(x)=="Presence")]
    })
    papa<-lapply(x@tesdat,function(x){
      y<-x[,"Presence"]
    })
    
  } else if(x@meta$replicatetype=="none"){
    
    outerloop<-1
    testa<-list(tester[,-which(colnames(tester)=="Presence")])
    papa<-list(tester[,"Presence"])
    
  }

  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="evaluation")

  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){

      # Make prediction
      pred=prd(x@fits[[i]][[j]],testa[[i]])
      
      #Feed with external threshold if available
      if(length(thres)==0){
        
        scores<-ceval(f=pred,
                      pa=papa[[i]],
                      tesdat=testa[[i]],
                      crit=crit)
        
      } else{
        
        scores<-ceval(f=pred,
                      pa=papa[[i]],
                      tesdat=testa[[i]],
                      tre=thres[which(names(thres)==names(x@fits[[i]])[j])],
                      crit=crit)
        
      } 
        
      lisa[[j]]<-scores

   }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa
  }
  names(lis)<-names(x@fits)
  out@performance<-lis
  
  return(out)
}

### =========================================================================
### define summary function for wsl.evaluation objects
### =========================================================================

sm=setMethod("summary",signature(object="wsl.evaluation"),definition=function(object){
  
  cat("\nMeta information: \n")
  df=data.frame(object@meta[c("author","date")],object@meta$wsl.fit[c("project","replicatetype","replicates")])  
  
  rownames(df)=""
  print(df)
  
  cat("\nThreshold: \n")
  df=as.data.frame(object@meta[c("cutoff")])  
  
  rownames(df)=""
  print(df)
  
  cat("\nMean skill: \n")
  
  mats=list()
  for(i in 1:length(object@performance)){
    mats[[i]]=do.call("cbind",object@performance[[i]])
  }
  mn=Reduce("+", mats) / length(mats)
  
  print(mn)
  
})

