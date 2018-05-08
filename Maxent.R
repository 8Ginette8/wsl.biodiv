# ###########################################################################
# Run Maxent

# $Date: 2014-06-18 V1
# $Date: 2018-03-14 V2
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
# Description: Flexible function to fit Maxent presence-only models.
# Required options:
# - feature types
# ...

# ###########################################################################
# Remove old Maxent Folders
library(dismo)
library(rJava)

# Run Maxent Model. Inputs:
# df.pres has to be a data.frame
# df.bg has to be a data.frame with the same columns as df.pres
# feat has to be a character indicating the features included

run.me<-function(dfra,type){
  
  if(type=="simple"){
    feat=c("linear","quadratic","hinge")
  } else if (type=="complex"){
    feat=c("linear","quadratic","hinge","product","threshold")
  }
 
  # Create directory for temporary MaxEnt data
  id <- as.numeric(commandArgs(trailingOnly=TRUE))
  me.temp.dir<-paste("Temp/tmp_Maxent",gsub("\\:","\\.",Sys.time()),sep="_")
  me.temp.dir<-paste(gsub(" ","_",me.temp.dir),id,sep="_")
  dir.create(me.temp.dir)
  
  # Translate features input into Maxent Command
  fiit.null<-c("linear","quadratic","hinge","product","threshold")
  fiit<-vector()
  for (i in fiit.null){
    if (is.element(i,feat)==T){
      fiit<-append(fiit,paste(i,"true",sep="="))
    }
    else{fiit<-append(fiit,paste(i,"false",sep="="))
    }
  }
  
  # Prepare data for MaxEnt model
  d.in<-dfra[,-which(colnames(dfra)=="Presence")]
  vec<-dfra[,"Presence"]
  
  #run maxent 
  me<-maxent(x=d.in,p=vec,args=fiit,path=me.temp.dir)    
  
  #Remove Temporary folder for Maxent
  unlink(me.temp.dir,recursive=T)
  
  return(me)
}

