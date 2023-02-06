### ==================================================================
### Loop over each species to extract the D2/R2 of each predictor
### ==================================================================
#' Predictive power assessments of predictors
#'
#' Evaluate the predictive power of single spatial predictors for fitting
#' spatial observations or continous values. Predictors may be classify in
#' classes to separate the ouputs, and data to explain (Y input) may be binary
#' (e.g. presences/absences), discrete (e.g. diversity count) or continuous.
#' Fix and random effects may be used with generalized linear model (GLM) and
#' simple linear model (LM). Therefore, the function currently implement
#' D2.adj::glm(), R2.adj::lm(), pseudo-R2::glmer() and pseudo-R2::lmer(). D2 is
#' calculated via the ecospat package and pseudo-R2 with the MuMIn package for mixed
#' models. Quadratic terms may also be used for more flexible fits, and a
#' parallel argument is available if processing the data is too time consuming.
#'
#' @param points Numeric. Object of class 'matrix' or 'data frame' with two columns named
#' "x" and "y", or 'SpatialPoints', or 'list' of 'SpatialPoints'. "points" should
#' be of same length as "val", "species" and "weight"
#' @param val Numeric. Object of class 'vector' if "points" is of class 'matrix' or 'data.frame',
#' otherwise object of class 'list'. Binary observations if "glmMODE="binomial",
#' discrete if "glmMODE="poisson", continious if "mlinear=TRUE". Must be of same
#' length as "points", "species" and "weight"
#' @param species Character. Object of class 'vector' specifying the ID of your observations.
#' No mandatory use except when "points" is a 'list' of 'SpatialPoints'. If used,
#' must be of same length as "points", "val" and "weight"
#' @param ras Object of class 'RasterLayer' or 'list' of 'RasterLayer'. Layers must be
#' of exact same resolution and extent i.e. exactly the same number of cells
#' @param rasCLASS Object of class 'vector' to associate IDs to the raster(s) if
#' "ras" is a 'list' of 'RasterLayer'
#' @param mlinear 'TRUE' run lm(), 'FALSE' run glm()
#' @param mixEffect Default is FALSE. If 'TRUE' run mLinear with random effect (i.e. mix LM or GLM)
#' @param mixCat Character. Object of class 'vector' if "points" is of class 'matrix' or 'data.frame',
#' otherwise object of class 'list'. Use to assign categories to your observations. Those
#' are used for potential random effects. Must be of same length as "points", "val" and "species"
#' @param glmMODE glm() link function if "mlinear=FALSE" (defaut::"binomial" or "poisson")
#' @param weight Numeric. Object of class 'vector' if "points" is of class 'matrix' or 'data.frame',
#' otherwise object of class 'list'. Used for weighting your observations when running glm()
#' or lm(). Must be of same length as "points", "val" and "species"
#' @param poly If TRUE, glm() or lm() uses a polynomial quadratic transformation on "ras"
#' @param polyBRUT If TRUE, force "poly" parameter in cases NAs are found in "ras"
#' @param parallel If TRUE, parallelisation is active
#' @param parINFOS Path. Create a txt file to write parallelisation infos if "parallel=TRUE".
#' @param cores Number of cores if "parallel=TRUE"
#' @param ... Arguments passed on to the lm() or glm() function
#'
#' @return 'wsl.varPPower' object with n slots corresponding to n "rasCLASS". In each
#' slots: matrix of nrow::c("species" + mean + standard deviation) & ncol::"ras". Where NAs
#' occur, models could not correctly converged.
#' @author Yohann Chauvier
#' @examples
#' ### Data Preparation
#' 
#' rm(list = setdiff(ls(), lsf.str()))
#' 
#' # Load environmental rasters and assign random raster Class
#' data(exrst)
#' rCLASS = c(rep("Pollution",4),rep("Climate",4),rep("LandCover",4))
#'
#' # Create a list of rasters out of the rasterBRICK
#' rasterL = unstack(rst)
#'
#' # Load my binary observations species data
#' data(var_select_XYtest)
#' 
#' # Create a category vector
#' mixV = sample(LETTERS[1:3], 35743, replace=TRUE)
#'
#' ### wsl.varPPower(): example with a data.frame & rasters not in classes
#' 
#' PPower.DF = wsl.varPPower(points=sp.DF[,c("x","y")],
#'                       val = sp.DF$myPA,
#'                       species = sp.DF$spCODES,
#'                       ras = rasterL,
#'                       rasCLASS = NULL,
#'                       mlinear = FALSE,
#'                       glmMODE = "binomial",
#'                       weight = sp.DF$myWEIGHT,
#'                       poly = TRUE,
#'                       polyBRUT = TRUE,
#'                       parallel = FALSE,
#'                       cores = NULL)
#' 
#' ### wsl.varPPower(): Same example including random effect
#' PPower.DF = wsl.varPPower(points = sp.DF[,c("x","y")],
#'                       val = sp.DF$myPA,
#'                       species = sp.DF$spCODES,
#'                       ras = rasterL,
#'                       rasCLASS = NULL,
#'                       mlinear = FALSE,
#'                       mixEffect = TRUE,
#'                       mixCat = mixV,
#'                       glmMODE = "binomial",
#'                       weight = sp.DF$myWEIGHT,
#'                       poly = TRUE,
#'                       polyBRUT = FALSE,
#'                       parallel = FALSE,
#'                       cores = NULL)
#'
#' ### wsl.varPPower(): example with a data.frame & rasters in classes
#' PPower.DF = wsl.varPPower(points=sp.DF[,c("x","y")],
#'                         val = sp.DF$myPA,
#'                         species = sp.DF$spCODES,
#'                         ras = rasterL,
#'                         rasCLASS = rCLASS,
#'                         mlinear = FALSE,
#'                         glmMODE = "binomial",
#'                         weight = sp.DF$myWEIGHT,
#'                         poly = TRUE,
#'                         polyBRUT = TRUE,
#'                         parallel = FALSE,
#'                         cores = NULL)
#' 
#' ### wsl.varPPower(): Same example including random effect
#' PPower.DF = wsl.varPPower(points=sp.DF[,c("x","y")],
#'                         val = sp.DF$myPA,
#'                         species = sp.DF$spCODES,
#'                         ras = rasterL,
#'                         rasCLASS = rCLASS,
#'                         mlinear = FALSE,
#'                         mixEffect = TRUE,
#'                         mixCat = mixV,
#'                         glmMODE = "binomial",
#'                         weight = sp.DF$myWEIGHT,
#'                         poly = TRUE,
#'                         polyBRUT = TRUE,
#'                         parallel = FALSE,
#'                         cores = NULL)
#' 
#' ### wsl.varPPower(): example with SpatialPointsDataFrame & rasters in classes
#'
#' # Checking length of every species input
#' c(length(mySP),length(myPA),length(spCODES),length(myWEIGHT))
#' 
#' # Generate new random effect classes
#' mixL = lapply(sapply(mySP,length),function(x) sample(LETTERS[1:3],x, replace=TRUE))
#'
#' # Running the function
#' PPower.SPDF = wsl.varPPower(points = mySP[[1]],
#'                           val = myPA[[1]],
#'                           species = NULL,
#'                           ras = rasterL,
#'                           rasCLASS = rCLASS,
#'                           mlinear = FALSE,
#'                           mixEffect=TRUE,
#'                           mixCat=mixL[[1]],
#'                           glmMODE = "binomial",
#'                           weight = myWEIGHT[[1]],
#'                           poly = TRUE,
#'                           polyBRUT = TRUE,
#'                           parallel = FALSE,
#'                           cores = NULL)
#' 
#' @export
wsl.varPPower=function(points,val,species=NULL,ras,rasCLASS=NULL,mlinear=FALSE,mixEffect=FALSE,
					   mixCat=NULL,glmMODE="binomial",weight=1,poly=FALSE,polyBRUT=FALSE,
					   parallel=FALSE,parINFOS=NULL,cores=detectCores()/2,...)
{
  ### ==================================================================
  ### Error handlings
  ### ==================================================================

  # Check mix effects parameters
  if (mixEffect&is.null(mixCat)){
  	stop("'mixCat' for mixEffect=TRUE was not specified...!")
  }

  # Checking class of input rasters
  # Checking if raster have the same CRS
  if (!(class(ras) %in% "list")) {
    cond1=class(ras)[1] %in% "RasterLayer"
    cond2=TRUE
  } else {
    r.class=unlist(lapply(ras,function(x) class(x)))
    cond1=all(r.class=="RasterLayer")
    r.crs=unlist(lapply(ras,function(x) proj4string(x)))
    cond2=length(unique(r.crs))==1
  }

  if (cond1==FALSE) {
    stop("Input raster(s) should be of class 'RasterLayer' !")
  }
  if (cond2==FALSE) {
    stop("Input raster(s) should have the same CRS !")
  }

  # Checking if input xy are correct
  cond3=!(class(points)[1] %in% "SpatialPoints")
  cond4=!(ncol(points)==2 || all(colnames(points)%in%c("x","y")))

  if (cond3 & cond4) {
    stop("Supplied points should be of class 'SpatialPoints' or a data.frame/matrix with two columns named x and y!")
  }

  # Checking if raster and points have the same CRS
  if (any(class(points) %in% c("matrix","data.frame"))) {
    if (!(compareCRS(ras[[1]],CRS("+init=epsg:4326")))) {
      warning("Make sure that input raster(s) and supplied XY points have the same CRS")
    }
  } else {
    cond5=c(class(points) %in% "list",class(ras) %in% "list")

    crsP=ifelse(cond5[1],unique(sapply(points,function(x) proj4string(x))),proj4string(points))
    crsR=ifelse(cond5[2],unique(sapply(ras,function(x) proj4string(x))),proj4string(ras))

    if(!(compareCRS(crsP,crsR))) {
      stop("Input raster(s) and supplied 'SpatialPoints' should have the same CRS")
    }
  }

  ### ==================================================================
  ### Conditions for 'species' and 'rasCLASS' parameters
  ### ==================================================================

  # Checking if the "species" parameter is not NULL
  cond6=is.null(species) & any(class(points) %in% c("matrix","data.frame"))
  cond7=is.null(species) & try(class(points@proj4string)[1],silent=T) %in% "CRS"

  if (cond6 | cond7) {
    spN=1
  } else {

    # Checking Length of 'points' and 'species' parameters
    spN=unique(species)
    if (length(spN)==1 & try(class(points@proj4string)[1],silent=T) %in% "CRS") {
      stop("'species' parameter should be NULL or have more than one species...")
    }

    LL=sapply(list(species,val,weight),function(x) length(x))

    if (any(class(points) %in% c("matrix","data.frame"))) {
      if (!(all(LL==nrow(points)))) {
        stop("XY points, 'species', 'val' & 'weight' should be of same length...")
      }
    } else {
      if (!(all(LL==length(points))) & is.null(species)) {
        stop("'species' parameter cannot be 'NULL' if 'list' of 'SpatialPoints' is used...")
      } else if (!(all(LL==length(points)))) {
        stop("SpatialPoints, 'species', 'val' & 'weight' should be of same length...")
      }
    }
  }

  # Checking if the "r.class" parameter is not NULL and assigning classes
  if (is.null(rasCLASS)) {
    VARclass="No class"
  } else {

    # Checking length of 'ras' and 'r.class' parameters
    VARclass=unique(rasCLASS)
    if (length(VARclass)==1) {
      stop("'rasCLASS' parameter should NULL or have more than one class...")
    }
    if (!(length(rasCLASS)==length(ras))) {
      stop("'rasCLASS' & 'ras' parameters should be of same length...")
    }
  }

  ### ==================================================================
  ### Applying models
  ### ==================================================================

  D1VAR=list()

  # Looping over the raster classes
  for (i in 1:length(VARclass))
  {
    cat("####### Processing class ","'",VARclass[i],"'",
        " (",i," out of ",length(VARclass)," Classes)...",sep="","\n")

    # Extract rasters from class i
    if (is.null(rasCLASS)) {
      VARn=ras
    } else {
      Q2=which(rasCLASS %in% VARclass[i])
      VARn=ras[Q2]
    }

    # Looping over the species

    if (parallel) {

      cat("Parallelisation is TRUE with",cores,"cores","\n")

      cl <- makeCluster(cores,outfile=parINFOS)
      registerDoParallel(cl)

    } else {

      cat("Parallelisation is FALSE","\n")

      registerDoSEQ()
    }


    D2SP=foreach (j=1:length(spN),.export="eco.adj.D2.glm",
      .packages=c("doParallel","foreach","raster","sp","lme4","MuMIn")) %dopar%
    {
      cat("Processing species ","'",unique(species)[j],"'",
          " (",j," out of ",length(spN)," Species)...",sep="","\n")

      # Keeping for each species the right info
      # If no species infos
      if (is.null(species)) {
        param=list(points,val,weight,mixCat)
      } else {
        Q1=which(species %in% spN[j])

        # If points are XY coords
        if (any(class(points) %in% c("matrix","data.frame"))) {
          param=list(points[Q1,],val[Q1],weight[Q1],mixCat[Q1])

          # If points are SpatialPoints
        } else {
          param=list(points[Q1][[1]],val[Q1][[1]],weight[Q1][[1]],mixCat[Q1][[1]])
        }
      }

      # Extracting for every XY the values of all rasters from the class
      if (class(VARn) %in% "list") {

        val.X=foreach (k=1:length(VARn)) %do%
        {
          cat("extract",names(VARn[[k]]),"\n")

          out=raster::extract(VARn[[k]],param[[1]],method="simple")
          return(out)
        }

        # Simple extract in case of just one raster
      } else {
        val.X=list(raster::extract(VARn[[1]],param[[1]],method="simple"))
      }

      # do.call and assign names to columns
      df.x=do.call("cbind",val.X)
      names(val.X)=sprintf("val.X%d",1:length(val.X))

      # Write formula for each raster val.X
      if (poly) {
        if (!polyBRUT)
        {
          form.m=lapply(names(val.X),function(x) {
          	as.formula(paste("val",paste("poly(",x,",2)"),sep="~"))})
        } else {
          form.m=lapply(names(val.X),function(x) {
          	as.formula(paste("val",paste(x,"+ I(",x,"^2)"),sep="~"))})
        }
      } else {
        form.m=lapply(names(val.X),function(x) {
        	as.formula(paste("val",paste(x),sep="~"))})
      }

      # Additonal preliminary changes if random effects 
      if (mixEffect) {
      	mix.form=lapply(form.m, function(x) {
      		as.formula(paste(paste(x[2],x[1],x[3],sep=""),"(1 | mixcat)",sep="+"))})
      	form.m=mix.form
      	df.mod=data.frame(val=param[[2]],w=param[[3]],mixcat=param[[4]],val.X)
      	df.mod$mixcat=as.factor(df.mod$mixcat)

      } else {
      	df.mod=data.frame(val=param[[2]],w=param[[3]],val.X)
      }

      # Apply formula GLM or LM
      if (mlinear)
      {
      	if (mixEffect) {

      		cat("Calculating pseudo-R2 adj.mixlm...",sep="","\n")
        	DR2=sapply(form.m,function(x) {
        		try.lmer=try(lmer(x,weights=w,data=df.mod),...)
        		if (class(try.lmer)%in%"try-error") {
        			return(NA)
        		} else {
        			return(r.squaredGLMM(try.lmer)[2])
        		}})*100

      	} else {

      		cat("Calculating R2 adj.lm...",sep="","\n")
        	DR2=sapply(form.m,function(x) {
        		summary(lm(x,weights=w,data=df.mod,...))$adj.r.squared})*100
      	}
  
      } else {

      	if (mixEffect) {

      		cat("Calculating pseudo-R2 adj.mixglm...",sep="","\n")
        	DR2=sapply(form.m,function(x) {
        		try.glmer=try(glmer(x,weights=w,family=glmMODE,data=df.mod,...))
        		if (class(try.glmer)%in%"try-error") {
        			return(NA)
        		} else {
        			return(r.squaredGLMM(try.glmer)[2])
        		}})*100

      	} else {

      		cat("Calculating D2 adj.glm...",sep="","\n")
        	DR2=sapply(form.m,function(x) {
        		eco.adj.D2.glm(glm(x,weights=w,family=glmMODE,data=df.mod,...))})*100
      	}
      }

      # Store
      return(DR2)
    }

    if (parallel) {stopCluster(cl)}

    # Aggregating list and assigning infos
    names(D2SP)=unique(species)
    result=do.call("rbind",D2SP)

    # Rename columns
    if (!(ncol(result)==1)) {
      colnames(result)=sapply(VARn,function(x) names(x))
    } else {
      colnames(result)=names(VARn)
    }

    # Do the summary
    if (nrow(result)==1) {
      meanSummary=NULL
      sdSummary=NULL
      quantilE=NULL
      result=rbind(result,result)
      row.names(result)=c("meanSummary","sdSummary")
    } else {
      # Creating a new row with the PPower summary
      meanSummary=apply(result,2,mean,na.rm=TRUE)
      quantilE=apply(result,2,quantile,na.rm=TRUE)
      sdSummary=apply(result,2,sd,na.rm=TRUE)
    }

    # Store species result for each assigned class of predictors
    result2=rbind(result,meanSummary,sdSummary,quantilE)
    D1VAR[[i]]=result2
  }
  names(D1VAR)=VARclass

  # Return the results
  return(D1VAR)
}