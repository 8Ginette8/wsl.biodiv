### =========================================================================
### preparation function
### =========================================================================
#' Check input data, collect meta information, take care of data subsetting.
#' Called by model fitting functions.
#'
#' Not to be called directly by the user
#' @author Philipp Brun, Yohann Chauvier
#' @export
preps=function(env=parent.frame(),call){
  
  env=as.list(env)
  
  ### ----------------
  ### check input data
  ### ----------------
  
  if(env$replicatetype%in%c("none","cv","block-cv","splitsample")==FALSE){stop("Non-existing replicatetype!")}
  
  if(env$replicatetype=="block-cv" && is.na(env$strata)){stop("Stratum vector needed for block crossvalidation!")}
  
  if(env$replicatetype=="block-cv" && length(unique(env$strata))!=env$reps){stop("Stratum vector has wrong number of levels!")}
  
  if(env$replicatetype%in%c("cv","block-cv","splitsample")==TRUE && (is.na("reps") ==TRUE || env$reps<2)){stop("Give reasonalbe number of replicates!")}
  
  if(env$replicatetype=="none" && env$reps>1){env$reps=1;warning("Replicate type is 'none' but multiple replicates were chosen - replicates are set to 1")}
  
  if(is.null(env$quadPoints)){
    
    if(env$replicatetype %in% "block-cv" && length(env$strata)!=length(env$pa)){
      stop("Length of strata should be equal to the number of P/A...")
    }

  } else if (class(env$quadPoints) %in% "ppp"){
    ref.Lppp=nrow(env$pres)
    if(env$replicatetype %in% "block-cv" && length(env$strata)!=ref.Lppp) {
      stop("wsl.ppmO: Length of strata should be equal to the number of presence points...")
    }

  } else if (class(env$quadPoints) %in% "wsl.quads") {
    ref.Lppp=nrow(env$pres)+nrow(env$quadPoints@coords)
    if(env$replicatetype %in% "block-cv" && length(env$strata)!=ref.Lppp) {
      stop("wsl.ppmGlaso: Length of strata should be equal to the number of pres/quadrature points...")
    }
  }

  #check for file path
  if(env$save && is.na(env$project)){stop("supply project in which data should be saved!")}
  
  if(env$save && !(env$project%in%list.dirs(env$path))){stop(paste("Project directory not existing in",env$path,
                                                                   "- please create manually"))} 
  
  ### -------------------
  ### generate wsl.fit obj
  ### -------------------
  
  out<-wsl.fit()
  
  # store function call
  out@call<-call
  call.t=as.character(call)
  
  ### -------------------
  ### add meta.info
  ### -------------------
  
  m.i=list()
  
  m.i$author=Sys.info()[["user"]]
  m.i$date=Sys.time()
  m.i$replicatetype=env$replicatetype
  m.i$replicates=env$reps
  m.i$taxon=env$taxon
  if (!(identical(grep("ppmO",call.t),integer(0))))
  {
    m.i$env_vars=paste(names(env$env_vars),collapse=", ")
  } else {
    m.i$env_vars=paste(colnames(env$env_vars),collapse=", ")
  }
  m.i$project=env$project
  m.i$model_tag=env$mod_tag
  
  # Add step info if exists
  if("step"%in%names(env)){
    m.i$step=env$step
  }

  # Add which_poly if exists
  if("which_poly"%in%names(env)){
    m.i$which_poly=env$which_poly
  }

  # Add infos on if predictor value = factor
  if (!(identical(grep("ppmG",call.t),integer(0)))){
    obs.fact=sapply(env$env_vars,class)%in%"factor"
    quad.fact=sapply(env$quadPoints@Qenv[,-1],class)%in%"factor"
    if (any((obs.fact+quad.fact)%in%1)){
      stop("Only 'env_vars' or 'quadPoints' contains predictor values of class 'factor'...")
    } else {
       m.i$which_factor=obs.fact
    }
  }
  
  out@meta=m.i
  
  ### ----------------------
  ### partition observations
  ### ----------------------
  
  # partition observations according to replicate type
  if (!(identical(grep("ppmO",call.t),integer(0)))) {
    
    obs.pres=env$pres
    names(obs.pres)=c("x","y")
    
    # We keep environmental values
    toRas=stack(lapply(env$env_vars,raster))
    toData=extract(toRas,obs.pres)
    dat=cbind(data.frame(obs.pres),toData)

  } else if (!(identical(grep("ppmG",call.t),integer(0)))) {
    
    obs.pres=env$pres
    names(obs.pres)=c("x","y")

    # We keep environmental values
    o.dat=cbind(data.frame(obs.pres),Presence=1,env$env_vars)
    q.dat=cbind(env$quadPoints@coords,env$quadPoints@Qenv)
    dat=rbind(o.dat,q.dat)

  } else {
    dat=cbind(data.frame(Presence=env$pa),env$env_vars)
  }
  
  obschoice<-list()
  testing<-list()
  XY<-list()
  if(env$replicatetype=="none"){
    obschoice[[1]]<-dat
    testing[[1]]<-dat
    
  } else if (env$replicatetype=="splitsample"){
    for (i in 1:env$reps){
      chc=sample(1:nrow(dat),size=round(nrow(dat)*0.7),replace=FALSE)
      obschoice[[i]]<-dat[chc,]
      testing[[i]]<-dat[-chc,]
      XY[[i]]<-env$xy[-chc,]
    }
    
  } else if (grepl("cv",env$replicatetype)){
    
    if(env$replicatetype=="cv"){
      unistr=sample(1:5,size=nrow(dat),replace=TRUE)      
    } else {
      unistr=env$strata   
    }
    
    for (i in 1:env$reps){
      obschoice[[i]]<-dat[which(unistr!=sort(unique(unistr))[i]),]
      testing[[i]]<-dat[which(unistr==sort(unique(unistr))[i]),]
      XY[[i]]<-env$xy[which(unistr==sort(unique(unistr))[i]),]
    }
  }
  
  # Add testing data to wsl.fit obj + coords if no PPPM
  out@tesdat=testing
  if (identical(grep("ppm",call.t),integer(0))){
    out@coords=XY
  } else {
    out@coords=list()
  }
  
  # return objects for model fitting
  return(list(wslfi=out,train=obschoice))
}