### =========================================================================
### BOYCE INDEX
### =========================================================================
#' Boyce index model evaluation sub-function II
#'
#' Not to be called directly by the user
#'
#' @author 'boycei' from the 'ecospat' package
#' @export
boycei <- function(interval, obs, fit) {
  
  fit.bin <- fit
  obs.bin <- obs
  fit.bin[fit[] >= interval[1] & fit[] <= interval[2]] <- "i"
  fit.bin[fit.bin != "i"] <- 0
  obs.bin[obs[] >= interval[1] & obs[] <= interval[2]] <- "i"
  obs.bin[obs.bin != "i"] <- 0
  
  pi <- length(which(obs.bin == "i"))/length(obs)
  ei <- length(which(fit.bin == "i"))/length(fit.bin)
  fi <- pi/ei
  
  return(fi)
}