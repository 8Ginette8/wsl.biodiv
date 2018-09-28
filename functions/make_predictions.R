# ###########################################################################
# Evaluation function and class definition
#
# $Date: 2018-05-17
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
#
#############################################################################
### =========================================================================
### Set Class
### =========================================================================

wsl.prediction<-setClass("wsl.prediction",slots=c(meta="list", # Meta information
                                                  thres="numeric", # supply external threshold
                                                  predictions="list")) # conserve function call

### =========================================================================
### define peval function
### =========================================================================
# 
wsl.predict<-function(x,predat=data.frame(),thres=numeric()){
  
  ### ------------------------
  ### Check input and prepare 
  ### ------------------------ 
  
  if(nrow(predat)==0){
    stop("Prediction data missing!")
  } 

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
    }
  } 
  
  
  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------
  
  out<-preva.meta(type="prediction")
  
  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){
      
      pred=prd(x@fits[[i]][[j]],predat)
      
      # Convert to binary predictions if thresholds were supplied
      if(length(thres)>0){
        the.tre=thres[which(names(thres)==names(x@fits[[i]])[j])]
        pred[pred>=the.tre]<-1
        pred[pred<the.tre]<-0
      }
     
      lisa[[j]]<-pred
      
    }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa
  }
  names(lis)<-names(x@fits)
  out@predictions<-lis
  
  return(out)
}
  
