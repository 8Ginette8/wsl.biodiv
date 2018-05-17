# ###########################################################################
# Second order functions for model evaluation
#
# $Date: 2018-05-08

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 

# ###########################################################################

### =========================================================================
### evaluate core function (peval)
### =========================================================================

ceval<-function(f,pa,tesdat,crit,tre=numeric()){
    
    # If there are any presences in the evaluation data
    if(any(pa==1)){

      z<-prediction(f,pa)
      # AUC
      auc<-performance(z,measure="auc")@y.values[[1]]
      rmse=performance(z,measure="rmse")@y.values[[1]]

      # optimum Threshold for conversion into binary outputs
      prbos<-seq(from=0,to=1,length.out=length(pa))
      zz<-performance(z,measure="sens",x.measure="spec",prbe=100000)
      zzz<-performance(z,measure="fpr",x.measure="fnr",prbe=100000)
      all.tss=zz@x.values[[1]]+zz@y.values[[1]]-1
      all.ppv<-performance(z,measure="ppv",prbe=100000)
      all.acc<-performance(z,measure="acc",prbe=100000)
      pn=zz@x.values[[1]]*zzz@x.values[[1]]
      py=zz@y.values[[1]]*zzz@y.values[[1]]

      all.kappa=(all.acc@y.values[[1]]-py*pn)/(1-py*pn)


      if(crit=="max"){
        ppv=max(all.ppv@y.values[[1]],na.rm=T)
        tss=max(all.tss)
        acc=max(all.acc@y.values[[1]],na.rm=T)
        kappa=max(all.kappa)
      } else {
        
        if(crit=="pp=op"){
          
          thr=quantile(f,probs=1-mean(pa))
          wi=which.min(abs(thr-z@cutoffs[[1]]))
        
        }else if(length(tre)!=0){
          
          wi=which.min(abs(tre-z@cutoffs[[1]]))
        }
        
        ppv=all.ppv@y.values[[1]][wi]
        tss=all.tss[wi]
        acc=all.acc@y.values[[1]][wi]
        kappa=all.kappa[wi]

      }

    # Return evaluation metrics
    weg=c(auc=auc,rmse=rmse,ppv=ppv,tss=tss,acc=acc,kappa=kappa)  
  
    return(weg)
    }
}

### =========================================================================
### prediction-evaluation meta.info function
### =========================================================================

preva.meta=function(env=parent.frame(),type=character()){
  
  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------
  
  m.i=list()
  m.i$author=Sys.info()[["user"]]
  m.i$date=Sys.time()
  m.i$wsl.fit=env$x@meta
  
  # Generate pevaluate object
  if(type=="evaluation"){
    out<-wsl.evaluation()
    m.i$cutoff=env$crit
  } else {
    out<-wsl.prediction()    
  }
  
  #add Meta info
  out@meta<-m.i
  
  return(out)
  
}

### =========================================================================
### prediction-evaluation prediction function
### =========================================================================

prd=function(mod,tst){
  
  # Generate probabilistic precitions
  if("MaxEnt"%in%class(mod)){
    
    pred<-predict(mod,tst)
    
  } else if("glm"%in%class(mod)){
    
    pred<-predict(mod,
                  newdata=tst,
                  type="response")
    
  } else if("gbm"%in%class(mod)){
    
    pred<-predict(mod,
                  newdata=tst,
                  n.trees=mod$n.trees,
                  type="response")
    
  } else if("randomForest"%in%class(mod)){
    
    pred<-predict(mod,
                  newdata=tst,
                  type="prob")[,2]
  }
  
  # Convert to numeric
  pred<-as.numeric(pred)
  
  return(pred)
  
}


