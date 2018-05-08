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

wsl.evaluation<-setClass("wsl.evaluation",slots=c(author="character", # Who created the object
                                           date="Date", # When was it created
                                           fitting.data="character", # summary of pfit.2 obj
                                           evaluation.type="character", # how should be evaluated
                                           evaluation.matrix="matrix", # what should be evaluated
                                           thres="numeric", # supply external threshold
                                           performance="list")) # performance objects for models

### =========================================================================
### define peval function
### =========================================================================
# 
wsl.evaluate<-function(x,em=NA,tester=data.frame(),thres=numeric()){
  
  ### ------------------------
  ### Check input and prepare Evaluation data
  ### ------------------------ 
  
  if(x@replicatetype=="none" && nrow(tester)==0){
    stop("External testing data must be supplied for replicatetype 'none'")
  } else if(x@replicatetype%in%c("crossvalidate","splitsample")) {
    outerloop<-length(x@tesdat)
    testa<-lapply(x@tesdat,function(x){
      y<-x[,-which(colnames(x)=="Presence")]
    })
    papa<-lapply(x@tesdat,function(x){
      y<-x[,"Presence"]
    })
    
  } else if(replicatetype=="none"){
    
    if(!(x@taxon%in%colnames(tester))){stop(paste("Taxon names do not match! Testing taxon should be",x@taxon))}
    
    outerloop<-1
    testa<-list(tester[,-which(colnames(tester)==x@taxon)])
    papa<-list(tester[,x@taxon])
    
  }

  ### ------------------------
  ### generate peval and add easy info
  ### ------------------------

  # Generate pevaluate object
  out<-wsl.evaluation()
  
  #add Meta info
  out@date<-Sys.Date()
  out@author<-basename(Sys.getenv("USERPROFILE"))
  out@fitting.data<-summary(x)
  
  # add evaluation type
  out@evaluation.type<-ifelse(x@replicatetype=="none","extrapolation",x@replicatetype)
  
  #define models to evaluate (default all tested models)
  if(length(em==1)){em<-x@model.matrix
      }else if(all(dim(em)==c(2,4)) && all(em%in%c(0,1))) {
        colnames(em)<-colnames(mm)
        rownames(em)<-colnames(mm)
      }else{stop("Evaluation matrix is ill-defined!")}
  
  out@evaluation.matrix<-em
  
  ### ------------------------
  ### evaluating tresholds
  ### ------------------------
  
  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)

  if(length(thres)>0){
    if(sum(em,na.rm=T)!=length(thres)){stop("wrong number of thresholds supplied!")}
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
      
      if(em[i,j]==1){
        
        print(paste("Evaluating",rownames(em)[i],colnames(em)[j],"..."))

        lisbet<-list()
        # Loop over replicates
        for(k in 1:x@replicates){

          # Generate probabilistic precitions
          if(colnames(em)[j]=="Maxent"){
            
            pred<-predict(x@fits[[i]]$Maxent[[k]],testa[[k]])
            
          } else if(colnames(em)[j]%in%c("GAM","GLM")){
            
            pred<-predict(x@fits[[i]][[which(names(x@fits[[i]])==colnames(em)[j])]][[k]],
                          newdata=testa[[k]],
                          type="response")
          
          } else if(colnames(em)[j]=="BRT"){
              
            pred<-predict(x@fits[[i]]$BRT[[k]],
                          newdata=testa[[k]],
                          n.trees=3500,
                          type="response")
          }
          
          # Convert to numeric
          pred<-as.numeric(pred)
          
          #Feed with external threshold if available
          if(length(thres)==0){
            
            scores<-ceval(pred,papa[[k]],testa[[k]],x@taxon)
            
          } else{
            
            modnam<-paste0(rownames(em)[j],colnames(em)[i])
            scores<-ceval(pred,papa[[k]],testa[[k]],x@taxon,tre=thres[which(names(thres)==modnam)])
            
          } 
          
          lisbet[[k]]<-scores

        }
        lisa[[j]]<-lisbet
      } else {
       lisa[[j]]<-list()
      }
      
    }
    names(lisa)<-colnames(em)
    lis[[i]]<-lisa
  }
  names(lis)<-rownames(em)
  out@performance<-lis
  
  return(out)
}
  
