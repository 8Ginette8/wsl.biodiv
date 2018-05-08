# ###########################################################################
# Run Random Forests

# $Date: 2014-08-15 V1
# $Date: 2018-03-14 V2
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL

# ###########################################################################

run.brt<-function(dfra,type){
  
  # Make weight vector
  wi=which(dfra[,"Presence"]==1)
  wt=rep(1,nrow(dfra))
  wt[wi]<-round((nrow(dfra)-length(wi))/length(wi))
  
  # Run models
  if(type=="simple"){
    waud=gbm(data=dfra,
             formula=Presence ~ . ,
             distribution = "bernoulli",
             interaction.depth = 1,
             shrinkage=.01,
             n.trees = 3500,
             weights=wt,
             verbose=F)
    
  } else if(type=="complex"){
    waud=gbm(data=dfra,
             formula=Presence ~ . ,
             distribution = "bernoulli",
             interaction.depth = 10,
             shrinkage=.003,
             n.trees = 3500,
             weights=wt,
             verbose=F)
  }

  return(waud)
 
}

