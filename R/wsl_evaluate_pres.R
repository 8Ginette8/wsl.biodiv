#' Evaluate presence-only models
#'
#' Assess model skill metrics for presence-only models in a wsl.fit object. Currently
#' Boyce index is evaluated. Furthermore, the threshold applied is returned.
#'
#' @param x A wsl.ppm fit object
#' @param tester Optional. A data.frame with testing data. Only mandatory if replicatetype='none'
#' was chosen when models were fitted. Otherwise, used when evaluation against external dataset is needed.
#' Must be a data.frame with as columns in order: "x", "y", "Presence" ('0' and '1'), "CV" (numeric: chosen
#' cv-folds; if replicatetype='none' -> only '1') and associated environmental values (same ones as
#' for fitted models; i.e. same columns order and names). Note that categorical predictor values must
#' be of class factor. NB: Here, model evaluation will only be initiated for the new testing data. 
#' @param thres Vector of the same length as the number of reps in model fit object
#' @param pres_abs Logical. If TRUE, evaluation metrics of presence-only models is applied to a wsl.fit
#' object of presence-absence models
#' @param env_vars Same spatial layers used in the fitted model or an object of class 'data.frame' or 
#' 'matrix' defining a sample of the target layers by keeping same order for columns. If spatial layers
#' are used, object of class 'RasterStack' or RasterBrick'.
#' @param window Only when 'wsl.ppmO' used'. Same object of class 'owin' used for models (in developments)
#' @param log_trans Logical. Use only if pres_abs=FALSE. Should predictions be converted to logarithm
#' before evaluation? Prevent model evaluation errors.
#' @param speedup If env_vars is a 'RasterStack' or RasterBrick', should the boyce evaluation be faster?
#' If TRUE, the algortihm uses a sample of the environmental layers
#' @param env_samp If speedup=TRUE, how many environmental cells should be sampled with replacement?
#' Default is 50'000 samples. The sample may be smaller than requested because of NAs
#' @param bias_cov A numerical vector whose length equal the number of environmental layers/columns.
#' Only used when a bias covariate is implemented in the PPM calibration i.e. to fit species observations
#' with a potential spatial observer bias. Default is 1 for each variable, whereas designated bias
#' covariate(s) (i.e. 0) will be reset everywhere to zero in order to evaluate corrected predictions
#' @param ... Additional arguments supplied to ecospat.boyce function (package 'ecospat')
#' @return an object of class 'wsl.evaluation'. If in slot "performance" NA as thresholds are found,
#' it indicates a lack of convergence in the model tested, and so, a biased/invalid threshold
#' @author Yohann Chauvier, Philipp Brun
#' @examples
#' 
#' 
#' ### Load
#' 
#' data(AlpineConvention_lonlat)
#' data(exrst)
#' rst = rst[[1:6]]
#' data(xy_ppm)
#' mypoints = xy.ppm[,c("x","y")]
#' 
#' ### Define mask
#' 
#' maskR = mask(rst[[1]],shp.lonlat)
#' 
#' ### Run 'wsl.ppm.window' function
#' 
#' wind = wsl.ppm.window(mask = maskR,
#'                       val = 1,
#'                       owin = TRUE)
#' 
#' ### Define quadrature points for 'wsl.ppmGlasso'
#' 
#'    # Grid regular
#' quadG1 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = FALSE,
#'                         lasso = TRUE,
#'                         env_vars = rst)
#' 
#' ### Define quadrature points for 'wsl.ppmO'
#' 
#'    # Grid regular
#' quadO1 = wsl.quadrature(mask = maskR,
#'                         area.win = wind,
#'                         random = FALSE,
#'                         lasso = FALSE,
#'                         env_vars = NULL)
#' 
#' ### Define your environments
#' 
#'    # For 'wsl.ppmGlasso' (observations focus)
#' envG = raster::extract(rst,mypoints)
#' 
#'    # For 'wsl.ppmO' (study area focus)
#' envO = wsl.ppm.env(rst,maskR)
#' 
#' ### Modelling
#' 
#'    # 'wsl.ppmGlasso' (alpha = 0.5 => Elastic net, see package 'glmnet')
#'        # Complex PPPM lasso (poly = TRUE & lasso=TRUE)
#' 
#' lasso1 = wsl.ppmGlasso(pres = mypoints,
#'                        quadPoints = quadG1@Qenv,
#'                        asurface = raster::area(shp.lonlat)/1000,
#'                        env_vars = envG,
#'                        taxon = "species_eg1",
#'                        replicatetype = "cv",
#'                        reps = 5,
#'                        strata = NA,
#'                        save=FALSE,
#'                        project = "lasso_eg1",
#'                        path = NA,
#'                        type = "binomial",
#'                        poly = TRUE,
#'                        lasso = TRUE,
#'                        alpha = 0.5,
#'                        type.measure = "mse",
#'                        standardize = TRUE,
#'                        nfolds = 5,
#'                        nlambda = 100)
#' 
#'        # Simple PPPM non lasso (poly = FALSE & lasso=FALSE)
#' 
#' lasso2 = wsl.ppmGlasso(pres = mypoints,
#'                        quadPoints = quadG1@Qenv,
#'                        asurface = raster::area(shp.lonlat)/1000,
#'                        env_vars = envG,
#'                        taxon = "species_eg2",
#'                        replicatetype = "cv",
#'                        reps = 5,
#'                        strata = NA,
#'                        save = FALSE,
#'                        project = "lasso_eg2",
#'                        path = NA,
#'                        type = "binomial",
#'                        mask = maskR,
#'                        poly = FALSE,
#'                        lasso = FALSE)
#' 
#'    # 'wsl.ppmO'
#'        # Simple PPPM non lasso (same as above with poly = FALSE & lasso = FALSE)
#' 
#' form.Sppm = as.formula(paste("~",paste(names(envO),collapse="+")))
#' lasso3 = wsl.ppmO(pres = mypoints,
#'                   quadPoints = quadO1,
#'                   env_vars = envO,
#'                   window = wind,
#'                   taxon = "species_eg3",
#'                   replicatetype = "cv",
#'                   reps = 5,                      
#'                   strata = NA,
#'                   save = FALSE,
#'                   project = "lasso_eg3",
#'                   path = NA,
#'                   formula = form.Sppm)
#' 
#'        # Complex PPPM non lasso
#' 
#' form.Cppm = as.formula(paste("~",paste(paste0("poly(",names(envO),",2)"),collapse="+")))
#' lasso4 = wsl.ppmO(pres = mypoints,
#'                   quadPoints = quadO1,
#'                   env_vars = envO,
#'                   window = wind,
#'                   taxon = "species_eg3",
#'                   replicatetype = "cv",
#'                   reps = 5,                      
#'                   strata = NA,
#'                   save = FALSE,
#'                   project = "lasso_eg3",
#'                   path = NA,
#'                   formula = form.Cppm)
#' 
#' ### Evaluation
#' 
#'    # Example for 'wsl.ppmGlasso'
#' 
#' eval1 = wsl.evaluate.pres(x = lasso1,
#'                           tester = NULL,
#'                           env_vars = rst,
#'                           mask = maskR,
#'                           window = NULL,
#'                           thres = NULL)
#' 
#' eval2 = wsl.evaluate.pres(x = lasso2,
#'                           tester = NULL,
#'                           env_vars = rst,
#'                           mask = maskR,
#'                           window = NULL,
#'                           thres = 0.001)
#'    
#'    # Example for 'wsl.ppmO'
#' eval3 = wsl.evaluate.pres(x = lasso3,
#'                           tester = NULL,
#'                           env_vars = envO,
#'                           mask = maskR,
#'                           window = wind,
#'                           thres = NULL)
#' 
#' eval4 = wsl.evaluate.pres(x = lasso4,
#'                           tester = NULL,
#'                           env_vars = envO,
#'                           mask = maskR,
#'                           window = wind,
#'                           thres = NULL)
#' 
#' @export
wsl.evaluate.pres<-function(x,tester=data.frame(),env_vars,window=NULL,thres=numeric(),
  pres_abs=FALSE,log_trans=TRUE,speedup=FALSE,env_samp=5e3,bias_cov=NULL,...){

  ### ------------------------
  ### check if pres_abs parameter well set
  ### ------------------------

  if (grepl("wsl.ppm",x@call)[1] && pres_abs){
    stop("'pres_abs' must be FALSE. 'x' is a wsl.fit object of presence-only models...")
  } else if (!(grepl("wsl.ppm",x@call)[1]) && !(pres_abs)){
    stop("'pres_abs' must be TRUE. 'x' is a wsl.fit object of presence-absence models...")
  }

  ### ------------------------
  ### check predictors depending on type of presence-only/pseudo-absences models
  ### ------------------------

  if (!(class(env_vars)[1]%in%c("RasterStack","RasterBrick","list","data.frame","matrix"))){
    stop("Wrong class of predictors for evaluation of presences-only models...")
  }

  if ("ppm"%in%class(x@fits[[1]][[1]]) & !(class(env_vars)[1]%in%"list")) {
    stop("'env_vars' must be a 'list' of predictors of class 'im' for Point Process Model evaluation...")
  } else if (!("ppm"%in%class(x@fits[[1]][[1]])) & !(class(env_vars)[1]%in%c("RasterStack","RasterBrick","data.frame","matrix"))) {
     stop("'env_vars' must be of class 'RasterStack', 'RasterBrick' or 'data.frame' for presence-only models evaluation...")
  }
  
  ### ------------------------
  ### check tresholds
  ### ------------------------
  
  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if (length(thres)>0){
    if (length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
      }
  }

  ### ------------------------------------------- 
  ### Which levels are used by the models if factors present?
  ### -------------------------------------------

  Xmeta=x@meta
  if (!(is.null(Xmeta$which_factor)))
  {
    if (any(Xmeta$which_factor)){
      col.fact = x@tesdat[[1]][,4:ncol(x@tesdat[[1]])][Xmeta$which_factor]
      id.factor = lapply(col.fact,levels)
    } else {
      id.factor = NULL
    }
  }

  ### ------------------------
  ### Check testing data and prepare for evaluation
  ### ------------------------ 

  if (x@meta$replicatetype=="none" && nrow(tester)==0){
    stop("External testing data must be supplied for replicatetype 'none'")
  }

  if (nrow(tester)!=0) {

    # Test if enough cv-fold in tester for number of fit reps
    cv.ind=length(table(tester$CV))
    fit.ind=length(x@fits)
    if (cv.ind!=fit.ind){
      warning("CV in 'tester' different than number of fit reps...!")
    }

    # Replace all tesdat
    p.cv = which(colnames(tester)%in%c("CV","x","y"))
    if (pres_abs){
      x@tesdat = lapply(1:fit.ind,function(x) tester[tester$CV%in%x,-p.cv])
    } else {
      x@tesdat = lapply(1:fit.ind,function(x) tester[tester$CV%in%x,-p.cv[3]])
    }
  }

  if (x@meta$replicatetype%in%c("cv","block-cv","splitsample")) {

    # If pres_abs=TRUE, we go look for coordinates
    if (pres_abs){
        coord.id <- lapply(x@tesdat,function(w) w$Presence%in%1)
        if (!(length(x@coords)%in%0)){
          testa.m <- lapply(1:length(coord.id),function(w) x@coords[[w]][coord.id[[w]],])
        } else {
          testa.m=list()
        }
    } else {
        testa.m <- lapply(x@tesdat,function(x) x[x$Presence%in%1,c("x","y")])
    }
  }
  
  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="evaluation")

  # If env_vars is a matrix change it to data.frame
  if (class(env_vars)[1]%in%'matrix') {env_vars = as.data.frame(env_vars)}

  ### ------------------------
  ### If there is a raster bias_cov reset it to zero
  ### ------------------------

  if (!is.null(bias_cov)){

    if (pres_abs) {b.i=-1} else {b.i=-(1:3)}
    tesdat.temp = lapply(x@tesdat,function(x) x[,b.i])
      
    for (zz in 1:length(tesdat.temp))
    {
      if (nrow(tesdat.temp[[zz]])!=0){
        tesdat.temp[[zz]][,bias_cov%in%0] = 0
        x@tesdat[[zz]] = cbind(x@tesdat[[zz]][,-b.i],tesdat.temp[[zz]])
      }
    }
    
    if (class(env_vars)[1]%in%c("RasterStack","RasterBrick")) {

      if (!(nlayers(env_vars)%in%length(bias_cov))) {
        stop("Length of 'bias_cov' vector must be equal to the number of raster layers...")
      }

      c.indbias = which(bias_cov%in%0)
      for (z in 1:length(c.indbias)) {env_vars[[c.indbias[z]]][] = 0}

    } else {

      if (!(ncol(env_vars)%in%length(bias_cov))) {
       stop("Length of 'bias_cov' vector must be equal to the number of raster layers...")
      }

      env_vars[,which(bias_cov%in%0)] = 0

    }
  }

  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # loop over replicates
  for (i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for (j in 1:length(x@fits[[1]])){

      # Try to open only the observations
      pres.only = x@tesdat[[i]][x@tesdat[[i]]$Presence%in%1,]

      if (all(is.na(x@fits[[i]][[j]]))) {

        # In case where nothing where fitted with ppmlasso because model was bad overall
        warning("Skipping Boyce evaluation: no model fit found!")
        scores = c(boyce=NA,threshold=NA) 

      } else if (nrow(pres.only)%in%0) {

        # We skip, in case of no "1" in the testing fold
        warning("Skipping presence-only evaluation: no observations found for testing-fold = ",i,"!")
        scores = c(boyce=NA,threshold=NA) 

      } else {
        if (pres_abs) {

          if (!(class(env_vars)[1] %in% "data.frame")) {

            if (length(testa.m)%in%0 || speedup) {
              ras.nan = which(!is.na(env_vars[[1]][]))
              samp.index = sample(ras.nan,env_samp,replace=TRUE)
              env_vars = rbind(pres.only[,-1],env_vars[][samp.index,])
            }
              
          } else {
            env_vars = rbind(pres.only[,-1],env_vars)
          }

          pred=prd.pa(x@fits[[i]][[j]],env_vars)

          if (any(class(pred) %in% c("data.frame","numeric"))) {
            pred = data.frame(pred,Pres=0)
            pred[1:nrow(pres.only),"Pres"] = 1
            pred = na.omit(pred)
          }
         
        } else {

          pres.only$Presence=NULL

          # Find out if polynomial term were used
          txt.call = paste0(gsub(" ","",deparse(x@call)),collapse="")
          pol = grepl("poly=TRUE",txt.call)

          # Provide which variables were poly() + associated coefs
          Xcoefs=x@coefs[[i]]
          
          # Make prediction & obtain boyce output
          if (speedup) {
            pred=prd.pres(mod=x@fits[[i]][[j]],env_vars=env_vars,window=window,
              polly=pol,meta=Xmeta,coefs=Xcoefs,id.fact=id.factor,valid.pres=pres.only,env_samp=env_samp)
          } else {
            pred=prd.pres(mod=x@fits[[i]][[j]],env_vars=env_vars,window=window,
              polly=pol,meta=Xmeta,coefs=Xcoefs,id.fact=id.factor,valid.pres=pres.only)
          }

          # Make the log.trans works only for non-lasso
          if (log_trans && !(class(x@fits[[i]][[j]]) %in% "cv.glmnet")) {
            if (class(pred) %in% "data.frame") {
              pred[,"pred"]=log(pred[,"pred"])
            } else {
              pred[]=log(pred[])
            }
          }
        }
          
        # Apply Boyce
        if (class(pred)%in%"data.frame"){
          pres.prob=pred[pred$Pres%in%1,]
          boyce<-ecoBoyce(fit=pred[,"pred"],obs=pres.prob[,"pred"],PEplot=FALSE,...)
        } else {
          boyce<-ecoBoyce(fit=pred,obs=testa.m[[i]],PEplot=FALSE,...)
        }

        # Feed with external threshold if available
        if (length(thres)==0){

          F0=boyce$F.ratio[boyce$F.ratio!=0]
          HS0=boyce$HS[boyce$F.ratio!=0]
          thre=try(HS0[-(1:tail(which(F0<1),1))][1],silent=TRUE)
          if (class(thre)%in%"try-error"){thre=NA}
          scores<-c(boyce=boyce$Spearman.cor,threshold=as.numeric(thre))

        } else {
          scores<-c(boyce=boyce$Spearman.cor,threshold=as.numeric(thres))
        }
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