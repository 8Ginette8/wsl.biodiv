### =========================================================================
### define peval function
### =========================================================================
#' Make predictions
#'
#' Make predictions with all models from a wsl.fit object. If thresholds are supplied
#' binary predictions are made, otherwise continuous predictions are returned.
#'
#' @param x An object of class wsl.fit
#' @param thres Optional. Object of the same length as the number of replicates, or model
#' types if a mean is applied accross model types. Obtained with 'get_thres'
#' @param predat Same spatial layers used in the fitted model or an object of class 'data.frame' or 
#' 'matrix' defining a sample of the target layers by keeping same order for columns. If spatial layers,
#' when wsl.ppmGlasso' used, object of class 'RasterStack' or RasterBrick'. When 'wsl.ppmO' used, a list of
#' 'im' object. 'wsl.ppmO' must be used with bias_cov = NULL (in developments)
#' @param window Only when 'wsl.ppmO' used'. Same object of class 'owin' used for models
#' @param log_trans Logical. Should the predictions be converted to logarithm before converting to binary ?
#' Should be TRUE if log.trans was TRUE when using wsl.evaluate
#' @param raster Logical. Should the output be a list of rasters or matrix ?
#' @param bias_cov A numerical vector whose length equal the number of environmental layers/columns.
#' Only used when a bias covariate is implemented in the PPM calibration i.e. to fit species observations with
#' a potential spatial observer bias. Default is 1 for each variable, whereas designated bias covariate(s)
#' (i.e. 0) will be reset everywhere to zero in order to evaluate corrected predictions
#' @return Object of class wsl.prediction with slots for meta info, and model predictions
#' @author Yohann Chauvier, Philipp Brun
#' @examples
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
#'                        mask = maskR,
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
#' ### Thresholds
#' 
#' get_thres(eval1, mean = FALSE)
#' get_thres(eval1, mean = TRUE)
#' 
#' get_thres(eval2, mean = FALSE)
#' get_thres(eval2, mean = TRUE)
#' 
#' get_thres(eval3, mean = FALSE)
#' get_thres(eval3, mean = TRUE)
#' 
#' get_thres(eval4, mean = FALSE)
#' get_thres(eval4, mean = TRUE)
#' 
#' ### Predictions
#' 
#'    # 'wsl.ppmGlasso'
#' pred1 = wsl.predict.pres(x = lasso1,
#'                          predat = rst,
#'                          window = NULL,
#'                          mask = maskR,
#'                          thres = get_thres(eval1,mean=FALSE),
#'                          raster = TRUE)
#' 
#' par(mfrow=c(2,3))
#' sapply(1:5,function(x) plot(pred1@predictions[[x]][[1]]))
#' 
#' pred2 = wsl.predict.pres(x = lasso2,
#'                          predat = rst,
#'                          window = NULL,
#'                          mask = maskR,
#'                          thres = NULL,
#'                          raster = TRUE)
#' 
#' par(mfrow=c(2,3))
#' sapply(1:5,function(x) plot(pred2@predictions[[x]][[1]]))
#' 
#' pred3 = wsl.predict.pres(x = lasso1,
#'                          predat = rst,
#'                          window = NULL,
#'                          mask = maskR,
#'                          thres = get_thres(eval1,mean=TRUE),
#'                          raster = TRUE)
#' 
#' par(mfrow=c(2,3))
#' sapply(1:5,function(x) plot(pred3@predictions[[x]][[1]]))
#' 
#' pred4 = wsl.predict.pres(x = lasso2,
#'                          predat = rst,
#'                          window = NULL,
#'                          mask = maskR,
#'                          thres = NULL,
#'                          raster = FALSE)
#' 
#'    'wsl.ppmO'
#' 
#' pred5 = wsl.predict.pres(x = lasso3,
#'                          predat = envO,
#'                          window = wind,
#'                          mask = NULL,
#'                          thres = get_thres(eval3,mean=FALSE),
#'                          raster = TRUE)
#' 
#' par(mfrow=c(2,3))
#' sapply(1:5,function(x) plot(pred5@predictions[[x]][[1]]))
#' 
#' pred6 = wsl.predict.pres(x = lasso4,
#'                          predat = envO,
#'                          window = wind,
#'                          mask = NULL,
#'                          thres = get_thres(eval4,mean=FALSE),
#'                          raster = TRUE)
#' 
#' par(mfrow=c(2,3))
#' sapply(1:5,function(x) plot(pred6@predictions[[x]][[1]]))
#' 
#' @export
wsl.predict.pres<-function(x,thres=numeric(),predat=list(),
  window=NULL,log_trans=TRUE,raster=FALSE,bias_cov=NULL){
  
  ### ------------------------
  ### Check input and prepare 
  ### ------------------------ 

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      if(length(x@fits)!=length(thres))
      stop("Wrong number of thresholds supplied! Should be one threshold per CV or model type...")
    }
  } 
  
  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------
  
  out <- preva.meta(type="prediction")

  # Find out if polynomial term were used
  txt.call = paste0(gsub(" ","",deparse(x@call)),collapse="")
  pol = grepl("poly=TRUE",txt.call)

  # Provide which variables were poly()
  Xmeta = x@meta

  # Which levels are used by the models if factors present?
  if (any(Xmeta$which_factor)){
    col.fact = x@tesdat[[1]][,4:ncol(x@tesdat[[1]])][Xmeta$which_factor]
    id.factor = lapply(col.fact,levels)
  } else {
    id.factor = NULL
  }

  # If predat is a matrix change it to data.frame
  if (class(predat)%in%'matrix') {predat=as.data.frame(predat)}

  ### ------------------------
  ### If there is a raster bias_cov reset it to zero
  ### ------------------------

  if (!is.null(bias_cov)){
    
    if (class(predat)%in%c("RasterStack","RasterBrick")) {

      if (!(nlayers(predat)%in%length(bias_cov))) {
        stop("Length of 'bias_cov' vector must be equal to the number of raster layers...")
      }

      c.indbias=which(bias_cov%in%0)
      for (z in 1:length(c.indbias)) {predat[[c.indbias[z]]][]=0}

    } else {

      if (!(ncol(predat)%in%length(bias_cov))) {
       stop("Length of 'bias_cov' vector must be equal to the number of raster layers...")
      }

      predat[,which(bias_cov%in%0)] = 0
    }
  }

  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){

      # Associated coefs
      Xcoefs = x@coefs[[i]]

      # Predictions
      pred = try(prd.pres(mod=x@fits[[i]][[j]],env_vars=predat,window=window,
        polly=pol,meta=Xmeta,coefs=Xcoefs,id.fact=id.factor),silent=TRUE)
      
      if (class(pred)%in%"try-error") {pred=NA}

      # Convert with log or not ?
      if (log_trans && !(class(x@fits[[i]][[j]]) %in% "cv.glmnet")) {
         if (class(pred) %in% "data.frame") {
           pred=log(pred$pred)
         } else {
           pred[]=log(pred[])
         }
      } else {
         if (class(pred) %in% "data.frame") {
           pred=pred$pred
         }
      }
      
      # Convert to binary predictions if thresholds were supplied
      if (length(thres)>0){
        
        if (length(x@fits)==length(thres)){
          the.tre=thres[which(names(thres)%in%names(x@fits)[i])][[1]][j]
        } else if (length(x@fits[[1]])==length(thres)){
          the.tre=thres[which(names(thres)==names(x@fits[[i]])[j])]
        } else {
          the.tre=thres
        }
        
        if (is.na(the.tre)){
          pred[]=NA
        } else {
          if (the.tre>0) {
            pred[][which(pred[]<the.tre)]<-0
            pred[][which(pred[]>0)]<-1
          } else {
            pred[][which(pred[]>=the.tre)]<-1
            pred[][which(pred[]<1)]<-0
          }
        }
      }

      if (raster){lisa[[j]]<-pred} else {lisa[[j]]<-pred[]}
      
    }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa

  }
  names(lis)<-names(x@fits)
  out@predictions<-lis
  out@thres<-unlist(thres)
  
  return(out)
}
