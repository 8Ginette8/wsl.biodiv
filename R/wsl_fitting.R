#' Fit presence-absence models
#'
#' Flexibly fit various types of functions but let framework take care of resampling,
#' meta-info storage, and file saving. wsl.flex basically allows supplying any possible
#' model, however, there may be problems with prediction/evaluation for exotic functions.
#' Available model algorithms for fitting are: Generalized Linear Models (GLM), Generalized
#' additive models (GAM), Maximum Entropy (MaxEnt), Artificial Neural Networks (ANN),
#' Generalized Boosted regression Models/Boosted Regression Trees (GBM) and Random Forest (RF)
#'
#' @param pa Object of class 'vector' with presence/absence values
#' @param env_vars Object of class 'data.frame' with environmental predictors
#' @param taxon  Name of the taxon for which models are fitted
#' @param x Optional. Object of class 'wsl.pseudoabsences'. If used 'pa', 'env_vars' and 'taxon'
#' will be replaced by new set of values set in the object.
#' @param replicatetype (How) should replicates be generated? May be 'none', 'splitsample',
#' 'cv' or 'block-cv'
#' @param reps Number of replicates
#' @param strata A numeric vector of the same length as observations with integers separating
#' cross validation replicates. Only used when replicatetype='block-cv'
#' @param save  Should the model be saved in a structured way? (not implemented yet)
#' @param project Character indicating the name of the project within which the models are run
#' (later used to define saving directories)
#' @param path Where to save? (not implemented yet)
#' @param step For glms and gams only. Should the models be updated with the step function?
#' @param mod_tag Not in wsl.flex. Descriptive label for current model
#' @param mod_args List with elements of class 'multi.input' which specify models to be fitted
#' in wsl.flex
#' @param xy Optional. XY coordinates of observations. Used for post evaluation with
#' presence-only metric when pres_abs=TRUE in wsl.evaluate.pres(). Must be a data.frame()
#' or matrix().
#' @return Object of class wsl.fit including slots for meta info, testing data for
#' evaluation, and model objects
#' @author Philipp Brun, Yohann Chauvier
#' @name fitdoc
#' @examples
#'
#' # Take anguilla data set from dismo package
#' data("Anguilla_train")
#' vrs=c("SegSumT","USRainDays","USSlope")
#' env=Anguilla_train[,vrs]
#'
#' ### Check out wsl.glm
#' form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
#'
#' modi1=wsl.glm(pa=Anguilla_train$Angaus,
#'               env_vars = env,
#'               taxon="Angaus",
#'               replicatetype="cv",
#'               reps=5,
#'               project="prototest",
#'               mod_tag="test-glm",
#'               formula=form.glm,
#'               family="binomial",
#'               step=TRUE)
#'
#' # Try out custom summary function
#' summary(modi1)
#'
#' # Access glm object of first replicate
#' summary(modi1@fits$replicate_01$`test-glm`)
#'
#' # Evaluate the model
#' eval1=wsl.evaluate.pa(modi1)
#'
#' # Get evaluation summary
#' summary(eval1)
#'
#' ### Check out wsl.gam
#' form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
#'
#' # Try out wsl.glm funcion
#' modi2=wsl.gam(pa=Anguilla_train$Angaus,
#'               env_vars = env,
#'               taxon="Angaus",
#'               replicatetype="splitsample",
#'               reps=3,
#'               project="prototest",
#'               mod_tag="test-gam",
#'               formula=form.gam,
#'               family="binomial",
#'               step=FALSE)
#'
#' # Try out custom summary function
#' summary(modi2)
#'
#' # Access glm object of first replicate
#' summary(modi2@fits$replicate_01$`test-gam`)
#'
#' # Evaluate the model
#' eval2=wsl.evaluate.pa(modi2,crit="maxTSS")
#'
#' # Get evaluation summary
#' summary(eval2)
#'
#' ### Check out wsl.gbm
#' form.gbm=as.formula(Presence ~ .)
#'
#' # Try out wsl.glm funcion
#' modi3=wsl.gbm(pa=Anguilla_train$Angaus,
#'               env_vars = env,
#'               taxon="Angaus",
#'               replicatetype="none",
#'               reps=1,
#'               project="prototest",
#'               mod_tag="test-brt",
#'               formula= form.gbm,
#'               distribution = "bernoulli",
#'               interaction.depth = 1,
#'               shrinkage=.01,
#'               n.trees = 3500)
#'
#' # Try out custom summary function
#' summary(modi3)
#'
#' # Access glm object of first replicate
#' summary(modi3@fits$replicate_01$`test-brt`)
#'
#' # Prepare external testing data
#' tste=data.frame(Presence=Anguilla_train$Angaus,CV=1,env)
#'
#' # Evaluate the model
#' eval3=wsl.evaluate.pa(modi3,crit="maxTSS",tester=tste)
#'
#' # Get evaluation summary
#' summary(eval3)
#'
#' ### Check out wsl.maxent
#' feat=c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")
#'
#' # Try out wsl.glm funcion
#' modi4=wsl.maxent(pa=Anguilla_train$Angaus,
#'                  env_vars = env,
#'                  taxon="Angaus",
#'                  replicatetype="block-cv",
#'                  reps=3,
#'                  strata=sample(1:3,nrow(env),replace=TRUE),
#'                  project="prototest",
#'                  mod_tag="test-mxe",
#'                  args=feat)
#'
#' # Try out custom summary function
#' summary(modi4)
#'
#' # Access glm object of first replicate
#' summary(modi4@fits$replicate_01$`test-mxe`)
#'
#' # Define external threshold
#' thmxe=c(`test-mxe`=0.5)
#'
#' # Evaluate the model
#' eval4=wsl.evaluate.pa(modi4,crit="external",thres=thmxe)
#'
#' # Get evaluation summary
#' summary(eval4)
#'
#' # Get thresholds
#' thr.4=get_thres(eval4)
#'
#' ### Check out wsl.flex
#' form.glm.2=as.formula(paste("Presence~",paste(vrs,collapse="+")))
#'
#' modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE,weight=TRUE),
#'    multi("gbm",list(formula=form.gbm,
#'    distribution = "bernoulli",
#'    interaction.depth = 1,
#'    shrinkage=.01,
#'    n.trees = 3500),"gbm-simple"), 
#'    multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE,weight=TRUE),
#'    multi("maxent",list(args=feat),"mxe-simple"),
#'    multi("randomForest",list(formula=form.gbm,ntree=500,maxnodes=NULL),"waud1"),
#'    multi("glm",list(formula=form.glm.2,family="binomial"),"glm-lin",step=TRUE,weight=TRUE))
#'
#' # Try out wsl.glm funcion
#' modi5=wsl.flex(pa=Anguilla_train$Angaus,
#'                env_vars = env,
#'                taxon="Angaus",
#'                replicatetype="block-cv",
#'                reps=3,
#'                strata=sample(1:3,nrow(env),replace=TRUE),
#'                project="multitest",
#'                mod_args=modinp)
#'
#' # Try out custom summary function
#' summary(modi5)
#'
#' # Access glm object of first replicate
#' summary(modi5@fits$replicate_01$`glm-simple`)
#'
#' # Evaluate the model
#' eval5<-wsl.evaluate.pa(modi5,crit="pp=op")
#'
#' # Get evaluation summary
#' summary(eval5)
#'
#' # Get thresholds
#' thr.5=get_thres(eval5)
#' 
#' @rdname fitdoc
#' @export
wsl.glm<-function(x=numeric(),
                  pa=numeric(),
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  step=FALSE,
                  mod_tag="",
                  xy=NULL,
                  ...){

  # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  #loop over replicates
  fits=list()
  for(i in 1:reps){
    modi=list()
    modi[[1]]=glm(...,data=lis$train[[i]])
    if(step){
      modi[[1]]=stepAIC(modi[[1]],direction="both",trace=FALSE)
    }
    names(modi)=ifelse(mod_tag=="","glm",mod_tag)
    fits[[i]]<-modi

  }

  # Name the fits
  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}

### =========================================================================
### define wsl.gam function
### =========================================================================
#' @rdname fitdoc
#' @export
wsl.gam<-function(x=numeric(),
                  pa=numeric(),
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  step=FALSE,
                  mod_tag="",
                  xy=NULL,
                  ...){

   # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  #loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    modi[[1]]=gam(...,data=lis$train[[i]])
    if(step){
      modi[[1]]=step(modi[[1]],direction="both",trace=FALSE)
    }
    names(modi)=ifelse(mod_tag=="","gam",mod_tag)
    fits[[i]]<-modi
  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}

### =========================================================================
### define wsl.maxent function
### =========================================================================
#' @rdname fitdoc
#' @export
wsl.maxent<-function(x=numeric(),
                     pa=numeric(),
                     env_vars=data.frame(),
                     taxon=character(),
                     replicatetype=character(),
                     reps,
                     strata=NA,
                     save=FALSE,
                     project=NA,
                     path=NA,
                     mod_tag="",
                     xy=NULL,
                     ...){

   # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()

    # Create directory for temporary MaxEnt data
    hde(me.temp.dir<-paste("tmp_Maxent",print(as.numeric(Sys.time())*1000,digits=15),sep="_"))
    dir.create(me.temp.dir)

    d.in<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
    vec<-lis$train[[i]][,"Presence"]

    modi[[1]]=maxent(x=d.in,p=vec,...,path=me.temp.dir)

    names(modi)=ifelse(mod_tag=="","mxe",mod_tag)
    fits[[i]]<-modi

    # Remove Temporary folder for Maxent
    unlink(me.temp.dir,recursive=TRUE)
  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}


### =========================================================================
### define wsl.gbm function
### =========================================================================
#' @rdname fitdoc
#' @export
wsl.gbm<-function(x=numeric(),
                  pa=numeric(),
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  mod_tag="",
                  xy=NULL,
                  ...){

  # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    modi[[1]]=gbm::gbm(...,data=lis$train[[i]])

    names(modi)=ifelse(mod_tag=="","gbm",mod_tag)
    fits[[i]]<-modi
  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}


### =========================================================================
### define wsl.ann function
### =========================================================================
#' @rdname fitdoc
#' @export
wsl.ann<-function(x=numeric(),
                  pa=numeric(),
                  env_vars=data.frame(),
                  taxon=character(),
                  replicatetype=character(),
                  reps,
                  strata=NA,
                  save=FALSE,
                  project=NA,
                  path=NA,
                  mod_tag="",
                  xy=NULL,
                  ...){

  # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  #loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    modi[[1]]=neuralnet(...,data=lis$train[[i]])
    names(modi)=ifelse(mod_tag=="","ann",mod_tag)
    fits[[i]]<-modi
  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}

### =========================================================================
### define wsl.flex function
### =========================================================================
#' @rdname fitdoc
#' @export
wsl.flex<-function(x=numeric(),
                   pa=numeric(),
                   env_vars=data.frame(),
                   taxon=character(),
                   replicatetype=character(),
                   reps,
                   strata=NA,
                   save=FALSE,
                   project=NA,
                   path=NA,
                   mod_args=list(),
                   xy=NULL){

  # Check supplied model types
  for(i in 1:length(mod_args)){
    if(!(mod_args[[i]]@mod%in%c("glm","gam","gbm","ann","maxent","randomForest"))){
      warning(paste(mod_args[[i]]@mod,"not in focal model functions. You might run in to problems when evaluating/predicting..."))
    }
  }

  # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    # loop over models
    for(j in 1:length(mod_args)){

      if(mod_args[[j]]@mod=="maxent"){

        # Create directory for temporary MaxEnt data
        hde(mod_args[[j]]@args$me.temp.dir<-paste("tmp_Maxent",
                                                  print(as.numeric(Sys.time())*1000,digits=15),
                                                  sep="_"))
        dir.create(mod_args[[j]]@args$me.temp.dir)

        mod_args[[j]]@args$x<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
        mod_args[[j]]@args$p<-lis$train[[i]][,"Presence"]

        hde(modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args))

        #Remove Temporary folder for Maxent
        unlink(mod_args[[j]]@args$me.temp.dir,recursive=T)

      } else {

        mod_args[[j]]@args$data=lis$train[[i]]

        # Make weight vector
        wi=which(mod_args[[j]]@args$data$Presence==1)
        wt=rep(1,nrow(mod_args[[j]]@args$data))
        wt[wi]<-round((nrow(mod_args[[j]]@args$data)-length(wi))/length(wi))

        if(mod_args[[j]]@weight){
          mod_args[[j]]@args$weights=wt
        }

        if(mod_args[[j]]@mod=="randomForest"){
          mod_args[[j]]@args$data$Presence=as.factor(mod_args[[j]]@args$data$Presence)
        }

        modi[[j]]=do.call(mod_args[[j]]@mod,mod_args[[j]]@args)
      }

      names(modi)[j]=ifelse(mod_args[[j]]@tag=="",paste0("model_",j),mod_args[[j]]@tag)

      if(mod_args[[j]]@step){
        modi[[j]]=stepAIC(modi[[j]],direction="both",trace=FALSE)
      }
    }

    fits[[i]]=modi

  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}
