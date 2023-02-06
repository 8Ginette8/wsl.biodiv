### =========================================================================
### define multi function
### =========================================================================
#' Define settings for function that should be supplied to wsl.flex
#'
#' Create a multi.input object that efficiently stores model specifications
#'
#' @param mod A character with the name of the function to be called. E.g. "gam"
#' @param args A list with arguments to be passed to the function specified in 'mod'
#' @param tag Character with name for model set-up
#' @param step Should step function be applied to update model
#' @return Object of class 'multi'
#' @author Philipp Brun
#' @examples
#'
#' ### Preliminary
#' 
#' data("Anguilla_train")
#' vrs=c("SegSumT","USRainDays","USSlope")
#' 
#' form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
#' form.gbm=as.formula(Presence ~ .)
#' form.glm.2=as.formula(paste("Presence~",paste(vrs,collapse="+")))
#' feat=c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")
#' form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
#' 
#' ### Multi examples
#'
#' multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE)
#' multi("gbm",list(formula=form.gbm,
#'                  distribution = "bernoulli",
#'                  interaction.depth = 1,
#'                  shrinkage=.01,
#'                  n.trees = 3500),"gbm-simple")
#' multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE)
#' multi("maxent",list(args=feat),"mxe-simple")
#' multi("randomForest",list(formula=form.gbm,ntree=500,maxnodes=NULL),"waud1")
#' multi("glm",list(formula=form.glm.2,family="binomial"),"glm-lin",step=TRUE)
#' 
#' @export
multi=function(mod,args,tag="",step=FALSE,weight=FALSE){

  out=multi.input()
  out@tag=tag
  out@step=step
  out@args=args
  out@mod=mod
  out@weight=weight

  return(out)
}
