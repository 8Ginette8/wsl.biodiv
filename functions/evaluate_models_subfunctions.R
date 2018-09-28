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


      if(crit=="maxTSS"){
        thr=z@cutoffs[[1]][which.max(all.tss)]

      } else if(crit=="pp=op"){
        thr=quantile(f,probs=1-mean(pa))
        
      }else if(length(tre)!=0){
        thr=tre
        
      }

    wi=which.min(abs(thr-z@cutoffs[[1]]))
    ppv=all.ppv@y.values[[1]][wi]
    tss=all.tss[wi]
    acc=all.acc@y.values[[1]][wi]
    kappa=all.kappa[wi]

    # Return evaluation metrics
    weg=c(auc=auc,rmse=rmse,ppv=ppv,tss=tss,acc=acc,kappa=kappa,threshold=as.numeric(thr))  
  
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
    
    pred<-df_or_rast(mod,tst)
    
  } else if("glm"%in%class(mod)){
    
    pred<-df_or_rast(mod=mod,
                  nwdat=tst,
                  type="response")
    
  } else if("gbm"%in%class(mod)){
    
    pred<-df_or_rast(mod,
                  nwdat=tst,
                  n.trees=mod$n.trees,
                  type="response")
    
  } else if("randomForest"%in%class(mod)){
    
    pred<-df_or_rast(mod,
                  nwdat=tst,
                  type="prob")
  }
  
  # Convert to numeric
  if(class(tst)=="data.frame"){
    pred<-as.numeric(pred)
  }
  
  return(pred)
  
}

### =========================================================================
### predict df or raster
### =========================================================================

df_or_rast=function(mod,nwdat,...){

  if("randomForest"%in%class(mod)){
    index=2
  }else{
    index=1
  }

  
  if(class(nwdat)%in%c("RasterStack","RasterBrick","RasterLayer")){
    
    add.arg=list(...)
    beginCluster(5)
    cl=getCluster()
    clusterExport(cl,list=list("nwdat","mod","add.arg","index"),envir=environment())
    out=clusterR(nwdat,predict,args=c(list(model=mod,index=index),add.arg))
    
    returnCluster()
    # out=raster::predict(nwdat,mod,...)
    endCluster()
    
  } else if(class(nwdat)=="data.frame"){
    out=predict(mod,newdata=nwdat,...)
    
    if("randomForest"%in%class(mod)){
      out=out[,2]
    }
  }

  return(out)
  
}
