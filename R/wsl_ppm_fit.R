#' Fit Poisson Point Process Models (PPPM)
#'
#' PPPM (Poisson Point Process Models) is a modelling approach genrally design to use
#' observation-only data aka point occurences. It allows the user to model species observations
#' intensity per unit area (i.e. the density of presence points over a spatial grid) as a log
#' linear function of predictors. Quadrature points (or background points) are necessary in the
#' model, and may be generated with the 'wsl.quadrature' function. Those points apply a spatial
#' scaling proportional to the study area and estimate the maximised model log likelihood (see
#' Renner 2013, Renner et al. 2015). Two version of the function are available. 'wsl.ppmO' applies
#' a simple point process model using the 'ppm' function from the package "spatstat" i.e. without
#' any further validation possible (only for calibration and prediction use). 'wsl.ppmGlasso'
#' applies a similar point process model but by implementing choices of regularisation and
#' variable selection following the package "glmnet". 
#'
#' @param pres Object of class 'data.frame' or 'matrix'. Coordinates xy of Species observations.
#' @param env_vars If 'wsl.ppmGlasso' used, must be an object of class 'matrix' or 'data.frame'
#' with environmental predictor values. If 'wsl.ppmO' used, must be an object of class 'list' generated
#' with the wsl.ppm.env' function. Note that categorical predictor values must be of class factor for
#' both 'env_vars' and 'quadPoints@Qenv'
#' @param quadPoints If 'wsl.ppmGlasso' used, must a 'wsl.quads' object generated with 'wsl.quadrature'
#' function and 'lasso=TRUE'. If 'wsl.ppmO' used, must be a 'ppp' object generated with the
#' 'wsl.quadrature' function with 'lasso=FALSE'; here, may be one same object for each 'reps' or a
#' list of objects whose length equal to 'reps'
#' @param asurface Only when 'wsl.ppmGlasso' used.The surface of the study area in square kilometers
#' @param window Only when 'wsl.ppmO' used. An object of class 'owin' generated with the
#' 'wsl.ppm.window' function
#' @param poly If TRUE, PPPM fits a second order polynomial regression
#' @param which_poly Which predictors should be using polynomial terms? Use a binary vector that
#' specify which variables/predictors. Length of vector must be equal to the number of input variables.
#' '1' stands for poly=TRUE whereas '0' stands for poly=FALSE. Default is polynomial for all if poly=TRUE.
#' @param lasso If FALSE no regularisation is applied
#' @param taxon Name of the taxon for which models are fitted
#' @param replicatetype (How) should replicates be generated? may be 'none', 'splitsample',
#' 'cv' or 'block-cv'
#' @param reps Number of replicates
#' @param strata A numeric vector of the same length as observations + quadrature points with
#' integers assigning cross validation replicates. Only used when replicatetype='block-cv'.
#' Note: the vector must first integrate CV information for observation points.
#' @param save  Should the model be saved in a structured way? (not implemented yet)
#' @param project Character indicating the name of the project within which the models are run
#' (later used to define saving directories)
#' @param path Where to save? (not implemented yet)
#' @param mod_tag Descriptive label for current model
#' @param formula Equation of the desired fit
#' @param penalty.glmnet If 'lasso=TRUE', a binary vector that specify which variables/predictors
#' used to model should be shrinked. Length of vector must be equal to the number of input
#' variables. '1' stands for shrinkage whereas '0' stands for no shrinkage, i.e. the variable
#' will always be included in the model.
#' @param ... If 'wsl.ppmGlasso' used with lasso = TRUE, arguments passed on to the
#' cv.glmnet() function (package 'glmnet') use to apply a Lasso, Ridge or Elastic Net
#' regularisation. To notice that the package's argument 'penalty.factor' is not needed
#' here. If used, the parameter 'penalty.glmnet' must instead be filled for each 'env_vars'.
#' If lasso = FALSE, arguments passed on to the glm("poisson") function. If 'wsl.ppmO'
#' used, arguments passed on to the ppm() function (package 'spatstat').
#' @return Object of class wsl.fit including slots for meta info, testing data for
#' evaluation, and model objects
#' @author Yohann Chauvier, Philipp Brun
#' @name fitppm
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
#'                        quadPoints = quadG1$Edata,
#'                        asurface = raster::area(shp.lonlat)/1000,
#'                        env_vars = envG,
#'                        taxon = "species_eg1",
#'                        replicatetype = "cv",
#'                        reps = 5,
#'                        strata = NA,
#'                        save=FALSE,
#'                        project = "lasso_eg1",
#'                        path = NA,
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
#'                        quadPoints = quadG1$Edata,
#'                        asurface = raster::area(shp.lonlat)/1000,
#'                        env_vars = envG,
#'                        taxon = "species_eg2",
#'                        replicatetype = "cv",
#'                        reps = 5,
#'                        strata = NA,
#'                        save = FALSE,
#'                        project = "lasso_eg2",
#'                        path = NA,
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
#' @rdname fitppm
#' @export
### ==================================================================
### PPM lasso fitting function
### ==================================================================


wsl.ppmGlasso<-function(pres=data.frame(),
                  env_vars=matrix(),
                  quadPoints=wsl.quads(),
                  asurface=numeric(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  mod_tag="",
                  poly=TRUE,
                  which_poly=NULL,
                  lasso=TRUE,
                  penalty.glmnet=NULL,
                  ...){

  # check and prepare data and output
  lis = preps(call=match.call())

  #loop over replicates
  fits = list()
  for (i in 1:reps) {

    # Extract 0 and 1 separately for each fold
    ppm.records=lis$train[[i]]
    ppm.pres=ppm.records[ppm.records$Presence%in%1,]
    ppm.quad=ppm.records[ppm.records$Presence%in%0,]

    # Original observations
    pres.env=as.data.frame(ppm.pres[,-c(1,2)])
    quad.env=as.data.frame(ppm.quad[,-c(1,2)])

    # Join quadrature and presences and apply ppm weight
    allD=rbind(pres.env,quad.env)
    allD$p.wt = rep(1.e-6,nrow(allD))
    allD$p.wt[allD$Presence==0] = asurface/nrow(quad.env)
    allD$Presence = allD$Presence/allD$p.wt
    allD = na.omit(allD)
    p.wt = allD$p.wt
    allD$p.wt = NULL

    modi = list()

    # Apply the lasso model or not ?
    if (lasso) {
      # Poly or not ?
      if (poly){
         
        pol.form=paste0("poly(",names(allD)[-1],",2)")

        if (!is.null(which_poly)) {
          if (length(names(allD)[-1])!=length(which_poly)){
            stop("Argument 'which_poly' should be of same length as number of predictors...")
          }
          pol.form[!(which_poly%in%1)]=names(allD)[-1][!(which_poly%in%1)]
        }
          
        formula = as.formula(paste("allD$Presence~",paste(paste0(pol.form),collapse="+")))

        # If factors, model.matrix needs to be adjusted
        if (any(sapply(allD,class)%in%"factor")){

          # Prepare contrasts param. to keep all factors in glmnet
          n.fact=names(allD)[sapply(allD,class)%in%"factor"]
          pre.contr=lapply(n.fact,function(x) allD[,x])
          post.contr=lapply(pre.contr, function(x) contrasts(x,contrasts=FALSE))
          names(post.contr)=n.fact

          # Apply the model.matrix with contrasts
          form = model.matrix(formula,data=allD,contrasts.arg=post.contr)

        } else {
          form = model.matrix(formula,data=allD)
        }

        # Store the polynomial coefs
        lis$wslfi@coefs[[i]] =
        lapply(allD[,-1][,grepl("poly",pol.form)],function(x) {
          popol = poly(x,2)
          return(attr(popol,"coefs"))
        })
      
      } else {
        form = as.matrix(allD)
      }

      if (is.null(penalty.glmnet)){

      	lassoG = try(cv.glmnet(form[,-1],allD$Presence,family="poisson",weight=p.wt,...),silent=TRUE)
      
      } else {

        if (length(penalty.glmnet)!=length(pol.form)) {
          stop("Arguments 'which_poly' and 'penalty.glmnet' should be of same length...")
        }

        # Set right penalty vectors
        which_factor=sapply(allD[,-1],class)%in%"factor"
        penalty = unlist(lapply(1:length(penalty.glmnet), function(x){
          
          if (!(grepl("poly",pol.form)[x])){
            if (which_factor[x]){
              p.glmnet=rep(penalty.glmnet[x],length(levels(allD[,-1][,which_factor])))
            } else {
               p.glmnet=penalty.glmnet[x]
            }
          
          } else {
            p.glmnet=rep(penalty.glmnet[x],2)}
            return(p.glmnet)
        }))

      	lassoG = try(cv.glmnet(form[,-1],allD$Presence,family="poisson",
          weight=p.wt,penalty.factor=penalty,...),silent=TRUE)
      }
      
      if (identical(class(lassoG),"try-error")) {modi[[1]] = NA} else {modi[[1]] = lassoG}
      names(modi) = ifelse(mod_tag=="","ppmlasso",mod_tag)
    
    } else {

      # Poly or not ?
      if (poly){

        pol.form=paste0("poly(",names(allD)[-1],",2)")

        if (!is.null(which_poly)) {
          if (length(names(allD)[-1])!=length(which_poly)){
            stop("Argument 'which_poly' should be of same length as number of predictors...")
          }
          pol.form[!(which_poly%in%1)]=names(allD)[-1][!(which_poly%in%1)]
        }
          
        form = as.formula(paste("allD$Presence~",paste(paste0(pol.form),collapse="+")))

        # Store the polynomial coefs
        lis$wslfi@coefs[[i]] =
        lapply(allD[,-1][,grepl("poly",pol.form)],function(x) {
          popol = poly(x,2)
          return(attr(popol,"coefs"))
        })

      } else {
        form = as.formula(paste("allD$Presence~",paste(paste0(names(allD)[-1]),collapse="+")))
      }

      modi[[1]] = glm(form,family=poisson(),weights=p.wt,data=allD,...)
      names(modi) = ifelse(mod_tag=="","ppmglm",mod_tag)
    }

    fits[[i]] = modi
  }
  
  # Name the fits
  names(fits) = paste0("replicate_",sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits = fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}


#' @rdname fitppm
#' @export
### ==================================================================
### PPM fitting function
### ==================================================================


wsl.ppmO<-function(pres=data.frame(),
              quadPoints=ppp(),
              env_vars=list(),
              window=owin(),
              taxon=character(),
              replicatetype=character(),
              reps,                      
              strata=NA,
              save=FALSE,
              project=NA,
              path=NA,
              formula=NULL,
              ...){

  # check and prepare data and output
  lis = preps(call=match.call())

  # Check quadPoints and reps
  if (!(class(quadPoints)%in%"ppp") & !(length(quadPoints)%in%reps)){
      stop("Error: Number 'quadPoints' object and number of 'reps' different...")
  }

  # Extract quadratures if a list
  if (class(quadPoints)%in%"ppp") {
    quadPoints=rep(list(quadPoints),reps)
  }
  
  #loop over replicates
  fits=list()
  for(i in 1:reps) {
    modi = list()
    Q = wsl.ppm.qscheme(lis$train[[i]],window,quadPoints[[i]])
    modi[[1]] = ppm(Q, trend=formula, covariates=env_vars,...)
    names(modi) = "ppm"
    fits[[i]] = modi
  }
  
  # Name the fits
  names(fits) = paste0("replicate_",sprintf("%02d",1:reps))
  
  # supply fitted objects
  lis$wslfi@fits = fits
  
  # Save
  #...
  
  return(lis$wslfi)
  
}
