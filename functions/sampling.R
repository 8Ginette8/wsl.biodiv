# ###########################################################################
# Function to sample for cross validation and pseudoabsences
#
# $Date: 2018-05-17
#
# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
#
#############################################################################
### =========================================================================
### cv blocks
### =========================================================================

make_blocks<-function(nstrat=4,df=data.frame(),nclusters=nstrat*5,npoints=NA){
  
  ### ------------------------
  ### check input data
  ### ------------------------
  
  if(nrow(df)==0 & is.na(npoints)){
    stop("Please supply number of points if no data.frame is supplied")
  }
  
  if(nrow(df)==0){
    
    ### ------------------------
    ### do ordinary sampling if no strata are supplied
    ### ------------------------
    
    out.strat=sample(rep(1:nstrat,ceiling(npoints/nstrat)),size=npoints)
    
  } else {
    
    ### ------------------------
    ### do stratification based on kmedoid clustering
    ### ------------------------
    
    # check for reasonable number of boxes
    if(nrow(df)<10*nclusters){
      stop("Too many boxes required!")
    }
    
    if(ncol(df)==1){
      
      clist=as.numeric(cut(df[,1],breaks=quantile(df[,1],probs=0:(nclusters)/(nclusters)),right=F))
      
      
    } else {
     
      # Scale input data
      scd=apply(df,2,scale)
      
      # do kmedoid clustering
      kmed=pam(scd,k=nclusters,metric="euclidean")
      
      # get clusters
      clist=kmed$clustering
      
    }
    
    # sort obtained clusters
    tbl=sort(table(clist),decreasing = T)
      
    ### ------------------------
    ### regularly assign clusters to strata
    ### ------------------------
    
    # prepare strata layers
    grps=rep(list(numeric()),nstrat)
    
    # for the clusters with many observations
    # distribute them regularly among strata but keep
    # last six clusters for estimating most 
    # regular distribution
    
    if(length(tbl)>6){
      
      for(i in 1:(length(tbl)-6)){
        
        # determine to which stratum the cluster should be
        # added
        fl=(floor((i-1)/nstrat))
        
        
        if(fl%%2==0){
          j=round(1+nstrat*((i-1)/nstrat-(floor((i-1)/nstrat))))
        } else {
          j=(nstrat+1)-round(1+nstrat*((i-1)/nstrat-(floor((i-1)/nstrat))))
        }
        
        # add cluster
        grps[[j]]=append(grps[[j]],tbl[i])
        
      }
    }
    
    # prepare for optimal distribution of last 6 clusters
    vlis=factor(1:nstrat,levels=1:nstrat)
    prs=rep(list(vlis),min(length(tbl),6))
    sstab=tbl[max(1,(length(tbl)-5)):length(tbl)]
    
    # Run brute-forcing gridSearch obtimization
    srch=gridSearch(levels=prs,fun=optme,nms=as.vector(sstab),grps=grps,tot=sum(tbl))
    
    # pull out results
    wi=as.numeric(as.character(srch$minlevels))
    
    # combine results with predistributed clusters
    for(i in 1:length(grps)){
      grps[[i]]=append(grps[[i]],sstab[wi==i])
    }
    
    # define vector with output strata
    out.strat=rep(NA,nrow(df))
    for(i in 1:length(grps)){
      out.strat[which(as.character(clist)%in%names(grps[[i]]))]=i
    }
  }
 
  # return result
  return(out.strat)
  
}

### =========================================================================
### optimization function for cluster distribution
### =========================================================================

optme=function(x,nms,grps,tot){
  
  # determine number of bservations in each groups from initial step
  grp=sapply(grps,"sum")
  
  # aggregate remaining observations by suggested cluster grouping
  x=as.numeric(as.character(x))
  agg.vals=aggregate(nms,by=list(x),FUN="sum")
  
  for(i in 1:length(grp)){
    
    if(i%in%agg.vals$Group.1){
      grp[i]=agg.vals$x[which(agg.vals$Group.1==i)]+grp[i]        
    }
    
  }
  
  # Calculate difference from equal distribution
  pen=(grp-tot/length(grp))^2

  return(sum(pen))
  
}

### =========================================================================
### sampling proportional to presence observation distribution
### =========================================================================

prop.sampling=function(points,nsamples=1000,res=1,adj=.2){
  
  ### ------------------------
  ### check input data
  ### ------------------------
  
  if(ncol(points)!=2 || !all(colnames(points)%in%c("x","y"))){
    
    stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
  }
  
  ### ------------------------
  ### Prepare point pattern object
  ### ------------------------
  
  # get x range
  xrng=range(points$x)
  
  # Determine order of magnitude of 1% of x range
  buff=0.01*abs(diff(xrng))
  oom=-10:10
  myoom=oom[which.min(abs(log10(buff)-oom))]
  
  # Extend x and y range by that magnitude to define observational window
  xrng[1]=floor(xrng[1]/10^myoom)*10^myoom
  xrng[2]=ceiling(xrng[2]/10^myoom)*10^myoom
  
  yrng=range(points$y)
  yrng[1]=floor(yrng[1]/100)*100
  yrng[2]=ceiling(yrng[2]/100)*100
  
  owi=owin(xrange=xrng,yrange=yrng)
  
  # Define Point Pattern object to calculate 
  myppp=ppp(x=points[,"x"],y=points[,"y"],window = owi)
  
  ### ------------------------
  ### Generate 'im' object with density info
  ### ------------------------
  
  lo=res*c(1000,round(1000*(yrng[2]-yrng[1])/(xrng[2]-xrng[1])))
  
  x=seq(xrng[1],xrng[2],length.out = lo[1])
  y=seq(yrng[1],yrng[2],length.out = lo[2])
  
  dens=density(myppp,adjust=adj,xy=list(x=x,y=y))
  
  ### ------------------------
  ### Draw locations proportional to point density
  ### ------------------------
  
  vls=as.vector(dens$v)
  
  # Replace NA's with zero probability
  if(any(is.na(vls))){
    vls[which(is.na(vls))]=0
  }
  
  csu=cumsum(vls)/(max(cumsum(vls)))
  
  draws=runif(nsamples)
  
  pts=sapply(draws,function(x,csu){
    diffe=csu-x
    wi=which(diffe>0)
    out=wi[which.min(diffe[wi])]
    return(out)
  },csu=csu)
  
  indi=arrayInd(pts,.dim=rev(lo))
  pt.out=cbind(x[indi[,2]],y[indi[,1]])
  colnames(pt.out)=colnames(points)
  
  return(pt.out)
  
}
