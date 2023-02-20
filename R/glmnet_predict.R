### =========================================================================
### SDM predict function for glmnet outputs
### =========================================================================

#' SDM predict function for glmnet outputs
#'
#' This function can make predictions based on any type of glmnet models, i.e.
#' even complex model calibrated with polynomial terms. Additionally, the predict
#' function avoids errors when pedictions are made based on predictors whose
#' values equal all one single value (very useful to apply bias covariate correction).
#'
#' @param mod A glmnet model output.
#' @param predat Same spatial layers used in the fitted model or an object of class 'data.frame'
#' or  matrix' defining a sample of the target layers by keeping same order for columns.
#' If spatial layers, when 'wsl.ppmGlasso' used, object of class 'RasterStack' or RasterBrick'.
#' @return RasterLayer or Vector.
#' @author Yohann Chauvier
#' @examples
#' 
#' @export
#' 
wsl.predictGlasso=function(mod,predat)
{
    
  # Predicting
  if("cv.glmnet"%in%class(mod)){

   # mask.env is a raster or a data.frame ?
   if (!(class(mask.env) %in% "data.frame")) {
      # Rasters into matrix
      dpl=lapply(1:nlayers(mask.env),function(x) mask.env[[x]][])
      env=as.data.frame(do.call("cbind",dpl))
   } else {
      env=mask.env
   }

   # Predicting and back into rasters (for poly and not poly)
   # coef(mod,s="lambda.min")
   p=try(predict(mod,newx=as.matrix(env),s="lambda.min"),silent=TRUE)

   if (class(p)%in%"try-error"){

      # Create polynomial values for our rasters
      envbis=na.omit(env)
      envbisi=na.action(na.omit(env))
      names(envbis)=names(mask.env)

     # Manage in case of constant predictors
      bID=apply(envbis,2,function(x) length(unique(x))==1)
      if (all(!bID)) {

        form1=as.formula(paste("~",paste(paste0("poly(",names(envbis),",2)")[!bID],collapse="+")))
        form3=model.matrix(form1,data=envbis)
    
      } else {

        # For constant predictors
        QB=envbis[,bID]
        form0=as.formula(paste("~",paste(paste0("poly(",names(envbis),",2,raw=TRUE)")[bID],collapse="+")))
        form01=model.matrix(form0,data=QB)[,-1]

        # For non constant predictors
        form1=as.formula(paste("~",paste(paste0("poly(",names(envbis),",2)")[!bID],collapse="+")))
        form2=model.matrix(form1,data=envbis)
        form3=cbind(form2,form01)
      }   
    
      Mvierge=matrix(nrow=dim(env)[1],ncol=dim(env)[2]*2)
      if (!is.null(envbisi)) {Mvierge[-envbisi,]=form3[,-1]} else {Mvierge=form3[,-1]}

      # Predict poly
      p=predict(mod,newx=as.matrix(Mvierge),s="lambda.min")
   }

   if (!(class(mask.env) %in% "data.frame")) {
      pred=raster(mask.env[[1]])
      pred[]=p
   } else {
      pred=as.numeric(p)
   }
   return(pred)
  } else {
    print("Not a glmnet model...")
  }
}