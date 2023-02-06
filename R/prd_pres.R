### =========================================================================
### prediction-evaluation function for presence-only model
### =========================================================================
#' Correctly feed the predict functions depending on model type (ppm, glm, gbm, maxent...)
#'
#' Not to be called directly by the user
#' @author Yohann Chauvier
#' @export
prd.pres=function(mod,env_vars,window,polly,meta,coefs,id.fact=NULL,valid.pres=NULL,env_samp=NULL){
  
  # Generate probabilistic precisions
  if ("ppm"%in%class(mod)) {
    if (is.null(window)) {stop("'window' parameter is needed for ppm evaluation...")}
    if (!(class(env_vars) %in% "data.frame")) {
      if (is.null(env_samp)){
        pred = raster(predict.ppm(mod,window=window,covariates=env_vars,ngrid=dim(window)))
      
      } else {
          ras.nan = which(!is.na(env_vars[[1]][]))
          samp.index = sample(ras.nan,env_samp,replace=TRUE)
          all.env = rbind(valid.pres[,-c(1,2)],env_vars[samp.index])
          pred = predict.ppm(mod,window=NULL,covariates=all.env,ngrid=NULL)

        if (!is.null(valid.pres)){
          pred = data.frame(pred,Pres=0)
          pred[1:nrow(valid.pres),"Pres"] = 1
          pred = na.omit(pred)
        }
      }
      
    } else {
      all.env = rbind(valid.pres[,-c(1,2)],env_vars)
      pred = predict.ppm(mod,window=NULL,covariates=all.env,ngrid=NULL)
      
      if (!is.null(valid.pres)){
          pred = data.frame(pred,Pres=0)
          pred[1:nrow(valid.pres),"Pres"] = 1
          pred = na.omit(pred)
      }
    }
    
  } else {

    # Predicting
    if("cv.glmnet"%in%class(mod)){

      # env_vars is a raster or a data.frame ?
      if (!(class(env_vars) %in% "data.frame")) {

        if (is.null(env_samp))
        {
          # Rasters into matrix
          dpl = lapply(1:nlayers(env_vars),function(x) env_vars[[x]][])
          env_ras = as.data.frame(do.call("cbind",dpl))
          colnames(env_ras)=names(env_vars)

        } else {
          ras.nan = which(!is.na(env_vars[[1]][]))
          samp.index = sample(ras.nan,env_samp,replace=TRUE)
          env_ras = rbind(valid.pres[,-c(1,2)],env_vars[samp.index])
        }

      } else {
        env_ras = rbind(valid.pres[,-c(1,2)],env_vars)
        
      }

      # Make sure that values from raster/DF is converted to factor if any (reset right levels too)
      env_ras[is.na(env_ras)]=NA
      if (any(meta$which_factor)){

        for (i in length(which(meta$which_factor)))
        {
          # Convert the target column to factor
          fr = which(meta$which_factor)[i]
          env_ras[,fr] = as.factor(env_ras[,fr])

          # Extract missing levels from calibration data and add them in order
          env.lvl = levels(env_ras[,fr])
          here.lvl = id.fact[[i]]%in%env.lvl
          levels(env_ras[,fr]) =  c(id.fact[[i]][here.lvl],id.fact[[i]][!here.lvl])
          env_ras[,fr] = factor(env_ras[,fr],levels=id.fact[[i]])
        }
      }

      # Predicting and back into rasters (for polly and not polly)
      if (!polly) {

        p = predict(mod,newx=as.matrix(env_ras),s="lambda.min")

      } else {

        # Create polynomial values for our rasters
        envbis = na.omit(env_ras)
        envbisi = na.action(na.omit(env_ras))
        names(envbis) = names(env_ras)

        # If values to predict -> only NAs, we return an NA raster
        if (nrow(envbis)==0){
          if (!(class(env_vars) %in% "data.frame")){
            pred = raster(env_vars[[1]])
            pred[] = NA
          } else {
            pred = rep(NA,nrow(env_vars))
          }
          return(pred)
        }

        # General form
        pol.form = paste0("poly(",names(envbis),", 2)")

        # Modify it according to which_poly
        if (!is.null(meta$which_poly)) {
          pol.form[!(meta$which_poly%in%1)] = names(envbis)[!(meta$which_poly%in%1)]
        }

        # Add the polynomial coeficients
        pol.form[grepl("poly",pol.form)] = paste0(gsub(")",",",pol.form[grepl("poly",pol.form)]),
          "coefs=",sprintf("coefs[[%d]]",1:length(coefs)),")")

        # Write formula
        form.m = as.formula(paste("~",paste(paste0(pol.form),collapse="+")))

        # If factors, model.matrix needs to be adjusted
        if (any(meta$which_factor)){

          # Prepare contrasts param. to keep all factors in glmnet
          n.fact = names(envbis)[meta$which_factor]
          pre.contr = lapply(n.fact,function(x) envbis[,x])
          post.contr = lapply(pre.contr, function(x) contrasts(x,contrasts=FALSE))
          names(post.contr) = n.fact

          # Apply the model.matrix with contrasts
          form3 = model.matrix(form.m,data=envbis,contrasts.arg=post.contr)

        } else {

          form3 = model.matrix(form.m,data=envbis)
        }

        # Store
        Mvierge = matrix(nrow=dim(env_ras)[1],ncol=ncol(form3[,-1]))
        if (!is.null(envbisi)) {Mvierge[-envbisi,] = form3[,-1]} else {Mvierge = form3[,-1]}
        colnames(Mvierge)=colnames(form3)[-1]

        # Predict poly
        p = predict(mod,newx=as.matrix(Mvierge),s="lambda.min")
      }

      if (!(class(env_vars) %in% "data.frame") && is.null(env_samp))
      {
        pred = raster(env_vars[[1]])
        pred[] = p

      } else {

        pred = as.numeric(p)
        
        if (!is.null(valid.pres)){
          pred = data.frame(pred,Pres=0)
          pred[1:nrow(valid.pres),"Pres"] = 1
          pred = na.omit(pred)
        }
      }

    } else {

       if (polly){
          new.names = gsub("poly(.)","",names(mod[1]$coefficients))
          new.names = unique(gsub(",.*","",new.names)[-1])
       } else {
          new.names = names(mod[1]$coefficients)[-1]
       }

      if (!(is.null(env_samp)) && !(class(env_vars) %in% "data.frame")) {

          ras.nan = which(!is.na(env_vars[[1]][]))
          samp.index = sample(ras.nan,env_samp,replace=TRUE)
          env_vars = as.data.frame(env_vars[samp.index])
      }

      if (class(env_vars) %in% "data.frame") {
        
        env_ras = rbind(valid.pres[,-c(1,2)],env_vars)
        names(env_ras) = new.names
        pred = stats::predict(mod,env_ras,type="response")

        if (!is.null(valid.pres)){
          pred = data.frame(pred,Pres=0)
          pred[1:nrow(valid.pres),"Pres"] = 1
          pred = na.omit(pred)
        }

      } else {
        
        names(env_vars) = new.names
        pred = raster::predict(env_vars,mod,type="response")
      }
    }
  }
  return(pred) 
}
