### =========================================================================
### pull out mean thresholds from evaluation
### =========================================================================

#' Get threshold
#'
#' Extracts thresholds from wsl.evaluation objects and names them so they can be fed
#' to the wsl.predict function. At the moment only averages over replicates can be obtained.
#'
#' @param x An object of class wsl.evaluation
#' @param mean Logical. If TRUE, a mean is applied to all thresholds
#' @return numeric() or 'vector' with thresholds
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
#' # Evaluate the model
#' eval1=wsl.evaluate.pa(modi1)
#'
#' # Get thresholds
#' get_thres(eval1, mean = FALSE)
#' get_thres(eval1, mean = TRUE)
#' 
#' @export
get_thres=function(x,mean) {

  if (class(x@performance[[1]]) %in% "numeric") {
    out=sapply(x@performance,function(x) x["threshold"])
    names(out)=gsub(".threshold","",names(out))

    if (mean) {
      out=mean(out,na.rm=TRUE)
    }
  
  } else {
    out=lapply(x@performance,function(x)
      sapply(x,function(y) y["threshold"]))
    new.n=gsub(".threshold","",names(out[[1]]))
    for (i in 1:length(out)) names(out[[i]])=new.n

    if (mean) {
      thress=do.call("rbind",out)
      out=colMeans(thress,na.rm=TRUE)

    }
  }
  return(out)
}
