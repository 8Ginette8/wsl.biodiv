# ###########################################################################
# Second order functions for model evaluation
#
# $Date: 2018-05-08

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 

# ###########################################################################

### =========================================================================
### evaluate core function (peval)
### =========================================================================

ceval<-function(f,pa,tesdat,tagg,tre=numeric()){
    
    # If there are any presences in the evaluation data
    if(any(pa==1)){
      
      z<-prediction(f,pa)
      # AUC
      auc<-performance(z,measure="auc")@y.values
      
      # optimum Threshold for conversion into binary outputs
      prbos<-seq(from=0,to=1,length.out=length(pa))
      zz<-performance(z,measure="tpr",x.measure="fpr",prbe=100000)
      zzz<-performance(z,measure="ppv",x.measure="npv",prbe=100000)
      maxPPV<-max(zzz@y.values[[1]],na.rm=T)
      meanPPV<-mean(zzz@y.values[[1]],na.rm=T)
      # decide whether to estimate threshold rom roc curve or to use
      # provided threshold
      
      if(length(tre)==0){
        roc<-cbind(z@cutoffs[[1]],zz@x.values[[1]],zz@y.values[[1]])[-1,] #,zzz@x.values[[1]],zzz@y.values[[1]]
        
        if(class(roc)=="matrix"){
          colnames(roc)<-c("Thres","FPR","TPR")
          thres.1<-roc[which.min(abs(roc[,3]-(1-roc[,2]))),1]
          thres.2<-quantile(f,probs=1-mean(pa))
          threse<-list(thres.1,thres.2)        
        } else {threse<-NA}
      } else{
        
        threse<-list(tre)
        print("..using supplied threshold")
        
      }
      
      # TSS
      out<-lapply(threse,function(thres){
        if(!(length(threse)==1 && is.na(threse))){
          bin<-rep(0,length(f))
          bin[f>thres]<-1
          differ<-pa-bin
          
          evil<-data.frame(diff=differ,truth=pa,pred=bin)
          sens<-evil[which(evil$truth==1),]
          spec<-evil[which(evil$truth==0),]
          trueprespred<-evil[which(evil$pred==1),]
          trueabspred<-evil[which(evil$pred==0),]
          
          sensitivity<-nrow(sens[which(sens$diff==0),])/nrow(sens)
          specificity<-nrow(spec[which(spec$diff==0),])/nrow(spec)
          
          tpp<-nrow(trueprespred[which(trueprespred$diff==0),])/nrow(trueprespred)
          tap<-nrow(trueabspred[which(trueabspred$diff==0),])/nrow(trueabspred)
          
          tss<- sensitivity + specificity - 1
          itss<- tpp + tap - 1
          
          predprev<-mean(bin)
          trueprev<-mean(pa)
          raus<-c(tss,thres,sensitivity,specificity,tpp,tap,itss,predprev,trueprev)
          names(raus)<-c("TSS","Threshold","Sensitivity","Specificity","PPV",
                         "NPV","FlippedTrueSkill","PredictedPrevalence","TruePrevalence")
          
        } else {raus<-"Model predictions failed"}
        
        return(raus)
      })
      names(out)<-c("Min_diff_ss","pp=op")
      
      
      weg<-unlist(list(AUC=auc,MaxPPV=maxPPV,MeanPPV=meanPPV,out))
      
      
    } else {
      
      weg<-rep(NA,21)
    }

  return(weg)
}
