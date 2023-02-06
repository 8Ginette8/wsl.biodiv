### =========================================================================
### prediction-evaluation function for all models
### =========================================================================
#' Correctly feed the predict functions depending on model type (glm, gbm, maxent...)
#'
#' Not to be called directly by the user
#' @author Philipp Brun
#' @export
prd.pa=function(mod,tst,clust=FALSE){

  # Generate probabilistic precitions
  if("MaxEnt"%in%class(mod)){

    pred<-df_or_rast(mod,tst,clust=clust)

  } else if("glm"%in%class(mod)){

    pred<-df_or_rast(mod=mod,
                     nwdat=tst,
                     type="response",
                     clust=clust)

  } else if("gbm"%in%class(mod)){

    pred<-df_or_rast(mod,
                     nwdat=tst,
                     n.trees=mod$n.trees,
                     type="response",
                     clust=clust)

  } else if("randomForest"%in%class(mod)){

    pred<-df_or_rast(mod,
                     nwdat=tst,
                     type="prob",
                     clust=clust)
  }

  # Convert to numeric
  if(class(tst)=="data.frame"){
    pred<-as.numeric(pred)
  }

  return(pred)

}
