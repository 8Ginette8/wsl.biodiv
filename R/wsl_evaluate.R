#' Evaluate presence-absence models
#'
#' Assess several model skill metrics for all models in a wsl.fit object. Currently
#' AUC, RMSE, TSS, PPV, Accuracy, and Cohen's Kappa are evaluated. Furthermore, the
#' threshold applied is returned.
#'
#' @param x A wsl.fit object
#' @param tester Optional. A data.frame with testing data. Only mandatory if replicatetype='none'
#' was chosen when models were fitted. Otherwise, used when evaluation against external dataset is needed.
#' Must be a data.frame with as columns in order : "Presence" ('0' and '1'), "CV" (numeric: chosen
#' cv-folds; if replicatetype='none' -> only '1') and associated environmental values (same ones as
#' for fitted models; i.e. same columns order and names). Note that categorical predictor values must
#' be of class factor. NB: Here, model evaluation will only be initiated for the new testing data. 
#' @param thres Vector of the same length as the number of reps in model fit object. For wsl.flex
#' model outputs, thresholds have to be labelled with the same names provided to models.
#' @param crit Which threshold criterion should be considered? Currently 'pp=op'
#' (predicted prevalence = observed prevalence), 'maxTSS' (threshold yielding maximum TSS),
#' and 'external' (thresholds manually supplied) are possible
#' @param prevalence_correction logical. Should imbalanced presence/absence data be upsampled to
#' prevalence 0.5 for model evaluation.
#' @param pres_only Logical. If TRUE, evaluation metrics of presence-absence models is applied
#' to a wsl.fit object of presence-only models i.e. generated with wsl.ppmGlasso()
#' @param window Only when 'wsl.ppmO' used'. Same object of class 'owin' used for models (in developments)
#' @param log_trans Logical. Use only if pres_only=TRUE. Should predictions be converted to logarithm
#' before evaluation? Prevent model evaluation errors.
#' @param bias_cov A numerical vector whose length equal the number of environmental layers/columns.
#' Only used when a bias covariate is implemented in the model calibration i.e. to fit species obs. with
#' a potential spatial observer bias. Default is 1 for each variable, whereas designated bias covariate(s)
#' (i.e. 0) will be reset everywhere to zero in order to evaluate corrected predictions
#' @return An object of class 'wsl.evaluation'
#' @author Philipp Brun, Yohann Chauvier
#' @examples
#' 
#' # Take anguilla data set from dismo package
#' data("Anguilla_train")
#' vrs=c("SegSumT","USRainDays","USSlope")
#' env=Anguilla_train[,vrs]
#'
#' ### Check out wsl.gam
#' form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
#'
#' # Try out wsl.gam funcion
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
#' # Access gam object of first replicate
#' summary(modi2@fits$replicate_01$`test-gam`)
#'
#' # Evaluate the model
#' eval2=wsl.evaluate.pa(modi2,crit="maxTSS")
#'
#' # Get evaluation summary
#' summary(eval2)
#' 
#' @export
wsl.evaluate.pa<-function(x,tester=data.frame(),window=NULL,thres=numeric(),
  crit="pp=op",prevalence_correction=FALSE,pres_only=FALSE,log_trans=TRUE,bias_cov=NULL){

  ### ------------------------
  ### check if pres_only parameter well set
  ### ------------------------

  if (!(grepl("wsl.ppm",x@call)[1]) && pres_only){
    stop("'pres_only' must be FALSE. 'x' is a wsl.fit object of presence-absence models...")
  } else if (grepl("wsl.ppm",x@call)[1] && !(pres_only)){
    stop("'pres_only' must be TRUE. 'x' is a wsl.fit object of presence-only models...")
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
    if (crit!="external"){
      warning("Assuming you want external tresholds to be used - setting crit='external'!")
      crit="external"
    }
  }

  if (!(crit%in%c("pp=op","maxTSS","external"))){
    stop("Invalid threshold criterion chosen!")
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

  if (x@meta$replicatetype=="none" && nrow(tester)==0) {
    stop("External testing data must be supplied for replicatetype 'none'")
  }

  if(prevalence_correction){
    x@tesdat = lapply(x@tesdat,function(y){
      tdpres = y[which(y$Presence==1),]
      tdabs = y[which(y$Presence==0),]
      if(nrow(tdabs)<nrow(tdpres)){
        tdabs = tdabs[sample(1:nrow(tdabs),nrow(tdpres),replace=TRUE),]
      } else if(nrow(tdpres)<nrow(tdabs)){
        tdpres = tdpres[sample(1:nrow(tdpres),nrow(tdabs),replace=TRUE),]
      }
      return(rbind(tdpres,tdabs))
    })
  }


  if (nrow(tester)!=0){

    # Test if enough cv-fold in tester for number of fit reps
    cv.ind=length(table(tester$CV))
    fit.ind=length(x@fits)
    if (cv.ind!=fit.ind){
      warning("CV in 'tester' different than number of fit reps...")
    }

    # Depending on number of reps, set new tesdat!
    p.cv = which(colnames(tester)%in%c("CV","x","y","Presence"))
    testa = lapply(1:fit.ind,function(x) tester[tester$CV%in%x,-p.cv])
    papa = lapply(1:fit.ind,function(x) tester[tester$CV%in%x,"Presence"])

  } else if (x@meta$replicatetype%in%c("cv","block-cv","splitsample")){

    # In case of pres_only=TRUE
    if (pres_only) {

      # Extract observations
      testa = lapply(x@tesdat,function(x) x[,-c(1:3)])
      papa = lapply(x@tesdat,function(x) x[,"Presence"])

    } else {
      testa = lapply(x@tesdat,function(x) x[,-which(colnames(x)=="Presence"),drop=FALSE])
      papa = lapply(x@tesdat,function(x) x[,"Presence"])
    }
  }

  # Bias_cov or not?
  if (!is.null(bias_cov)){
    
    if (ncol(testa[[1]])!=length(bias_cov)){
      stop("Length of 'bias_cov' vector must be equal to the number of environmental columns/layers...")
    }

    for (zz in 1:length(testa)){
      if (nrow(testa[[zz]])!=0) {testa[[zz]][,bias_cov%in%0] = 0}
    }
  }

  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="evaluation")

  ### -------------------------------------------
  ### Evaluate models
  ### -------------------------------------------

  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){

    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){

      # We skip, in case of only '1' or '0' in the testing fold
      if (length(table(papa[[i]]))%in%1|length(table(papa[[i]]))==0){

        warning("Skipping presence-absence evaluation: no '0' or '1' found for testing-fold = ",i,".NAs are returned!")
        scores = c(auc=NA,rmse=NA,ppv=NA,tss=NA,acc=NA,kappa=NA,threshold=NA) 

      } else if (any(table(papa[[i]])<2)){

        warning("Skipping presence-absence evaluation: too few '0' or '1' found for testing-fold = ",i,".NAs are returned!")
        scores = c(auc=NA,rmse=NA,ppv=NA,tss=NA,acc=NA,kappa=NA,threshold=NA) 
      
      } else {

        # Make prediction
        if (pres_only) {

          # Find out if polynomial term were used
          txt.call = paste0(gsub(" ","",deparse(x@call)),collapse="")
          pol = grepl("poly=TRUE",txt.call)

          # Provide which variables were poly() + associated coefs
           if (pol) {
            Xcoefs = x@coefs[[i]]
          }
          Xmeta = x@meta
        
          pred = prd.pres(mod=x@fits[[i]][[j]],env_vars=testa[[i]],valid.pres=NULL,
            window=window,polly=pol,meta=Xmeta,coefs=Xcoefs,id.fact=id.factor)

          if (log_trans && !(class(x@fits[[i]][[j]]) %in% "cv.glmnet")) {
            pred = log(pred)
          }
       
        } else {
          pred = prd.pa(x@fits[[i]][[j]],testa[[i]])
        }

        # Feed with external threshold if available
        if(length(thres)==0){

          scores = ceval(f=pred,
                        pa=papa[[i]],
                        crit=crit)

        } else {

          scores = ceval(f=pred,
                        pa=papa[[i]],
                        tre=thres[which(names(thres)==names(x@fits[[i]])[j])],
                        crit=crit)

        }
      }

      lisa[[j]]<-scores

    }
    names(lisa) = names(x@fits[[i]])
    lis[[i]]<-lisa
  }
  names(lis)<-names(x@fits)
  out@performance<-lis

  return(out)
}
