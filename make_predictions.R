# ###########################################################################
# Fit decade wise niches
#
# $Date: 2014-06-24
#
# Author: Philipp Brun, pgbr@aqua.dtu.dk
# National Institute of Aquatic Resources (DTU-Aqua), 
# Technical University of Denmark, Charlottenlund
# 
# Description: This file contains the definition of the generic class "predis"
# as well as the associated function "prez".
# Goal of the predis class is to allow for a flexible and convenient predicting
#
#############################################################################
### =========================================================================
### Set Class
### =========================================================================

# generate pfit class
predis<-setClass("predis",slots=c(author="character", # Who created the object
                                 date="Date", # When was it created
                                 fitting.data="character", # summary of pfit.2 obj                               
                                 prediction.matrix="matrix", # what should be predicted
                                 thres="numeric", # supply threshold to factorize predictions
                                 predictions="list")) # Prediction vectors

### =========================================================================
### define peval function
### =========================================================================
# 
prez<-function(x,pm=NA,predat=data.frame(),thres=numeric()){
  
  ### ------------------------
  ### Check input and prepare 
  ### ------------------------ 
  
  if(nrow(predat)==0){
    stop("Prediction data missing!")
  } 
  
  if(x@replicatetype%in%c("crossvalidate","splitsample")) {
    outerloop<-length(x@tesdat)
  } else if(x@replicatetype=="none"){
    outerloop<-1
  }

  ### ------------------------
  ### generate peval and add easy info
  ### ------------------------

  # generate predis object
  out<-predis()
  
  # add meta info
  out@date<-Sys.Date()
  out@author<-basename(Sys.getenv("USERPROFILE"))
  out@fitting.data<-summary(x)
  
  #define models to evaluate (default all tested models)
  if(length(pm==1)){pm<-x@model.matrix
      }else if(all(dim(pm)==c(2,4)) && all(pm%in%c(0,1))) {
        colnames(pm)<-colnames(mm)
        rownames(pm)<-colnames(mm)
      }else{stop("Prediction matrix is ill-defined!")}
  
  out@prediction.matrix<-pm
  
  ### ------------------------
  ### evaluating tresholds
  ### ------------------------
  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)

  if(length(thres)>0){
    if(sum(pm,na.rm=T)!=length(thres)){stop("wrong number of thresholds supplied!")}
  }
  
  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # Loop over complexity
  for(i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){
      
      if(pm[i,j]==1){
        
        print(paste("Predicting",rownames(pm)[i],colnames(pm)[j],"..."))

        lisbet<-list()
        # Loop over replicates
        for(k in 1:x@replicates){

          # Generate probabilistic precitions
          if(colnames(pm)[j]=="Maxent"){
            
            # add workaround to preserve NA slots
            cc=which(complete.cases(predat))
            pva<-predict(x@fits[[i]]$Maxent[[k]],predat[cc,])
            pred=rep(NA,nrow(predat))
            pred[cc]<-pva
            
          } else if(colnames(pm)[j]%in%c("GAM","GLM")){
            
            pred<-predict(x@fits[[i]][[which(names(x@fits[[i]])==colnames(pm)[j])]][[k]],
                          newdata=predat,
                          type="response")
          
          } else if(colnames(pm)[j]=="BRT"){
              
            pred<-predict(x@fits[[i]]$BRT[[k]],
                          newdata=predat,
                          n.trees=3500,
                          type="response")
          }
          
          # Convert to numeric
          pred<-as.numeric(pred)
          
          # Convert to binary predictions if thresholds were supplied
          if(length(thres)>0){
            pred[which(pred>=thres)]<-1
            pred[which(pred<thres)]<-0
          }
          
          lisbet[[k]]<-pred

        }
        lisa[[j]]<-lisbet
      } else {
       lisa[[j]]<-list()
      }
      
    }
    names(lisa)<-colnames(pm)
    lis[[i]]<-lisa
  }
  names(lis)<-rownames(pm)
  out@predictions<-lis
  
  return(out)
}
  
