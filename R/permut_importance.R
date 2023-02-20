### =========================================================================
### Assess predictor importance using random permutations
### =========================================================================

#' Assess predictor importance using random permutations
#'
#' This function is based on the principle of permutation importance. This metric is applicable
#' to any type of model. The basic idea is to consider a variable important if it has a positive
#' effect on the prediction accuracy (classification), or MSE (regression). The risk of using this
#' metric is a potential bias towards collinear predictive variables. It is important to note that
#' this function is fully compatible with singular SDM output of the package.
#'
#' @param model A model output. e.g. GLM, GAM, GBM, random forest, glmnet outputs...
#' @param data A data.frame of environmental values.
#' @param nperm How many random permutations should be applied per predictor?
#' @param type Response type. Type of desired response.
#' @param rescale Should the results be rescaled from a scale from 1 to 100?
#' @return Average model importance score of each predictor.
#' @author Patrice Descombes
#' @examples
#'
#' # Take anguilla data set from dismo package
#' data("Anguilla_train")
#' vrs = c("SegSumT","USRainDays","USSlope")
#'
#' # Apply a simple GLM and measure each predictors' contribution to the model
#' form.glm = as.formula(paste("Angaus~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
#' glm.calib = glm(form.glm,data=Anguilla_train,family="binomial")
#' imp.test = var.imp(model = glm.calib,
#'                    data = Anguilla_train[,vrs],
#'                    nperm = 10,
#'                    type = "response",
#'                    rescale = TRUE);imp.test
#' 
#' # Apply a glmnet model and measure each predictors' contribution to the model
#' 
#' @export
#' 
var.imp <- function (model, data, nperm, type, rescale) {
  if (type=="response") {
    if (class(model)[1] == "gbm") {
      ref <- predict(model, data, type="response", n.trees=model$n.trees)
    } else if (class(model)[1] == "cv.glmnet") {
      ref <- wsl.predictGlasso(model,data)
    } else {
      ref <- predict(model, data, type="response")
    }
  }
  if (type=="prob") {
    ref <- predict(model, data, type="prob")[,2]
  }
  VarImp <- vector()
  for (i in 1:ncol(data)) {
    print(names(data)[i])
    refi <- vector()
    for (j in 1:nperm) {
      cali <- data
      cali[,i] <- cali[sample(1:nrow(cali), nrow(cali)), i]
      if (type=="response") {
        if (class(model)[1] == "gbm") {
          refi <- c(refi, 1 - cor(ref, predict(model, cali, type="response", n.trees=model$n.trees)))
        } else if (class(model)[1] == "cv.glmnet") {
          refi <- c(refi, 1 - cor(ref, wsl.predictGlasso(model,cali)))
        } else {
          refi <- c(refi, 1 - cor(ref, predict(model, cali, type="response")))
        }
      }
      if (type=="prob") {
        refi <- c(refi, 1 - cor(ref, predict(model, cali, type="prob")[,2]))
      }
    }
    VarImp <- c(VarImp, round(mean(refi), 3))
  }
  names(VarImp) <- names(data)
  if (rescale==TRUE) {
    VarImp <- VarImp/sum(VarImp)*100
  }
  return(VarImp)
}