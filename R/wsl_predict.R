### =========================================================================
### define peval function
### =========================================================================
#' Make predictions
#'
#' Make predictions with all models from a wsl.fit object. If thresholds are supplied
#' binary predictions are made if convert=TRUE, otherwise continuous predictions and
#' separate description of thresholds are returned
#'
#' @param x An object of class wsl.fit
#' @param predat Data.frame or raster for which predictions should be made
#' @param thres Optional. Object of the same length as the number of replicates, or model
#' types if a mean is applied accross model types. Obtained with 'get_thres'
#' @param bias_cov A numerical vector whose length equal the number of environmental layers/columns.
#' Only used when a bias covariate is implemented in calibrations i.e. to fit species observations with
#' a potential spatial observer bias. Default is 1 for each variable, whereas designated bias covariate(s)
#' (i.e. 0) will be reset everywhere to zero in order to evaluate corrected predictions
#' @param clust Logical. If raster predictions are made, should the operation be run in parallel ?
#' @return Object of class wsl.prediction with slots for meta info, and model predictions
#' @author Philipp Brun, Yohann Chauvier
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
#' thr.4=get_thres(eval4, mean=FALSE)
#'
#' ### Check out wsl.flex
#' form.glm.2=as.formula(paste("Presence~",paste(vrs,collapse="+")))
#'
#' modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE,weight=TRUE), 
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
#' thr.5=get_thres(eval5, mean=FALE)
#'
#' ### Make some predictions
#' pred4=wsl.predict.pa(modi4,predat=env)
#' pred5=wsl.predict.pa(modi5,predat=env,thres=thr.5)
#' 
#' @export
wsl.predict.pa<-function(x,predat=data.frame(),thres=numeric(),bias_cov=NULL,clust=FALSE) {
  
  ### ------------------------
  ### Check input and prepare 
  ### ------------------------ 
  
  if(nrow(predat)==0){
    stop("Prediction data missing!")
  } 

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      if(length(x@fits)!=length(thres))
      stop("Wrong number of thresholds supplied! Should be one threshold per CV or model type...")
    }
  }


  ### ------------------------
  ### Bias covariates 
  ### ------------------------ 


  # Bias_cov or not?
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
  
  
  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------
  
  out<-preva.meta(type="prediction")
  
  ### ------------------------------------------- 
  ### Evaluate models
  ### -------------------------------------------
  
  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){
    
    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){
      
      pred=prd.pa(x@fits[[i]][[j]],predat,clust)
      
      # Convert to binary predictions if thresholds were supplied
      if(length(thres)>0){

        if (length(x@fits[[1]])==length(thres)){
          the.tre=thres[which(names(thres)==names(x@fits[[i]])[j])]
        } else if (length(x@fits)==length(thres)){
          the.tre=thres[which(names(thres)%in%names(x@fits)[i])][[1]][j]
        }
         
        if (the.tre>0) {
          pred[][which(pred[]<the.tre)]<-0
          pred[][which(pred[]>0)]<-1
        } else {
          pred[][which(pred[]>=the.tre)]<-1
          pred[][which(pred[]<1)]<-0
        }
      }
     
      lisa[[j]]<-pred
      
    }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa

  }
  names(lis)<-names(x@fits)
  out@predictions<-lis
  out@thres<-unlist(thres)
  
  return(out)
}
