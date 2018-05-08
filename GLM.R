# ###########################################################################
# Run GLMs

# $Date: 2014-08-19 V1
# $Date: 2018-03-14 V2
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL

# ###########################################################################

# Run P/A GAM. Inputs:
# dfra: input data frame containing all data
# type: what type of GAM should it be (currently simple and complex exist)

run.glm<-function(dfra,type){
  
  ### ------------------------
  ### Generate Formula
  ### ------------------------ 

  # Get predictive part
  preds<-predictors(dfra,"Presence",type=type,mod="glm")
  
  # Make weight vector
  wi=which(dfra[,"Presence"]==1)
  wt=rep(1,nrow(dfra))
  wt[wi]<-round((nrow(dfra)-length(wi))/length(wi))
  
  # Make presence/absence formula
  foermeli<-as.formula(paste0("Presence~",preds))
  
  ### ------------------------
  ### Run GLM
  ### -----------------------
 
  glmmel<-glm(foermeli,family=binomial(link="logit"),data=dfra,weights=wt)

   #Stepwise model selection
  print("stepwise selecting GLM formula...")
  glmmel<-stepAIC(glmmel,direction="both",trace=F)

  return(glmmel)  
}

