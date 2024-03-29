### =========================================================================
### Subsample based on two stratification critera
### =========================================================================
#' Subsample based on two stratification critera
#' @param spp Spatial data.frame
#' @param type either 'env.strat' or 'env.semi.strat'
#' @param sampsiz desired size of the subsample
#' Not to be called directly by the user
#' @references
#' Descombes, P., Chauvier, Y., Brun, P., Righetti, D., Wüest, R. O., Karger, D. N., ... &
#' Zimmermann, N. E. (2022). Strategies for sampling pseudo-absences for species distribution
#' models in complex mountainous terrain. bioRxiv, 2022-03.
#' @author Philipp Brun
#' @export
stratify=function(spp,type,sampsiz){

  # Identify number of points per stratum
  spp=spp[sample(1:nrow(spp)),]

  tb=table(spp$layer)
  stb=sort(tb)
  cms=cumsum(stb)

  if(type=="env.strat"){
    # Identify sample size by including rare strata

    tarsiz=optimize(function(x,sz,tb){
      vls=cbind(tb,x)
      sm=sum(apply(vls,1,min))
      return(abs(sz-sm))
    },tb=tb,sz=sampsiz,interval=c(sampsiz/length(tb),5*sampsiz/length(tb)))

    gregr=aggregate(as.data.frame(coordinates(spp)),
                    by=list(lay=spp$layer),function(x,y){
                      out=x[1:min(length(x),y)]
                      return(out)
                    },y=ceiling(tarsiz$minimum))

  } else {

    rat=sampsiz/sum(log(tb))

    tarsiz=optimize(function(x,sz,tb){
        vls=cbind(log(tb)*x,tb)
        sm=sum(apply(vls,1,min))
        return(abs(sz-sm))
    },tb=tb,sz=sampsiz,interval=c(rat,5*rat))

    gregr=aggregate(as.data.frame(coordinates(spp)),
                    by=list(lay=spp$layer),function(x,y){
                      out=x[1:min(ceiling(log(length(x))*y),length(x))]
                      return(out)
                    },y=ceiling(tarsiz$minimum))


  }

  smpsiz=sapply(gregr[,2],function(x){
    length(x)
  })

  # if there is only one sample per layer just take as is
  if(is.numeric(gregr[,2])){
    bgstrt2=SpatialPointsDataFrame(gregr[,2:3],
                                   data=data.frame(layer=gregr[,1]),
                                   proj4string = spp@proj4string)
  } else {
    # otherwise reformat
    bgstrt2=SpatialPointsDataFrame(cbind(do.call("c",gregr[,2]),do.call("c",gregr[,3])),
                                   data=data.frame(layer=rep(gregr$lay,smpsiz)),
                                   proj4string = spp@proj4string)
  }



  return(bgstrt2)
}
