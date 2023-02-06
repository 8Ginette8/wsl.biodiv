#' Pre-select set of predictors
#'
#' Rank spatial predictors "ras" in descending order following their mean predictive
#' power found with wsl.varPPower(), and apply descending successive correlations tests
#' to obtain a pre-selection of predictors. Additionally a Variance Inflation Factor
#' Test (VIF) may also be applied to predictors which succeeded the descending tests to
#' identify the ones that show high mutlicollinearity.
#'
#' @param PPower.object 'wsl.varPPower' object
#' @param ras Object of class 'RasterLayer' or 'list' of 'RasterLayer'. Same as in wsl.varPPower().
#' Layers must be of same resolution i.e. excatly the same number of cells
#' @param corTEST Numeric. Threshold of accepted correlation applied in the descending tests
#' @param vifTEST Logical. If TRUE, a VIF test is applied to the set of predictors that succeeded
#' descending correlation tests
#' @param corVIF Numeric. Threshold of accepted correlation within VIF if "vifTEST=TRUE"
#' @param ... Arguments passed on to the cor() function used for "corTEST"
#' @return Object of class list with three elements: 1) ranking of non correlated predictors,
#' 2) ranking of leftovers correlated predictors, 3) vif test results 
#' @author Yohann Chauvier
#' @examples
#' ### Data Preparation
#'
#' # Load environmental rasters and assign random raster Class
#' data(exrst)
#' rCLASS = c(rep("Pollution",4),rep("Climate",4),rep("LandCover",4))
#'
#' # Create a list of rasters out of the rasterBRICK
#' rasterL = unstack(rst)
#'
#' # Load my binary observations species data
#' data(var_select_XYtest.RData)
#'
#' ### wsl.varPPower(): example with a data.frame & rasters not in classes
#' 
#' PPower.DF = wsl.varPPower(points=sp.DF[,c("x","y")],
#'                           val = sp.DF$myPA,
#'                           species = sp.DF$spCODES,
#'                           ras = rasterL,
#'                           rasCLASS = NULL,
#'                           mlinear = FALSE,
#'                           glmMODE = "binomial",
#'                           weight = sp.DF$myWEIGHT,
#'                           poly = TRUE,
#'                           polyBRUT = TRUE,
#'                           parallel = FALSE,
#'                           cores = NULL)
#'
#'
#' ### wsl.varPPower(): example with SpatialPointsDataFrame & rasters in classes
#' 
#' PPower.DF = wsl.varPPower(points=sp.DF[,c("x","y")],
#'                           val = sp.DF$myPA,
#'                           species = sp.DF$spCODES,
#'                           ras = rasterL,
#'                           rasCLASS = NULL,
#'                           mlinear = FALSE,
#'                           glmMODE = "binomial",
#'                           weight = sp.DF$myWEIGHT,
#'                           poly = TRUE,
#'                           polyBRUT = TRUE,
#'                           parallel = FALSE,
#'                           cores = NULL)
#'
#' ### wsl.varPPower(): example with SpatialPointsDataFrame & rasters in classes
#'
#' # Checking length of every species input
#' c(length(mySP),length(myPA),length(spCODES),length(myWEIGHT))
#'
#' # Running the function
#' PPower.SPDF = wsl.varPPower(points = mySP,
#'                             val = myPA,
#'                             species = spCODES,
#'                             ras = rasterL,
#'                             rasCLASS = rCLASS,
#'                             mlinear = FALSE,
#'                             glmMODE = "binomial",
#'                             weight = myWEIGHT,
#'                             poly = TRUE,
#'                             polyBRUT = TRUE,
#'                             parallel = FALSE,
#'                             cores = NULL)
#'
#'
#'  ### wsl.varKeep() with a descending correlation of 0.7 with a VIF test at 0.7
#' 
#' # Example with no raster classes
#' KeepVar = wsl.varKeep(PPower.object = PPower.DF,
#'                       ras = rasterL,
#'                       corTEST = 0.7,
#'                       vifTEST = TRUE,
#'                       corVIF = 0.7,
#'                       use="complete.obs")
#'
#' #Example with raster classes
#' KeepVar = wsl.varKeep(PPower.object = PPower.SPDF,
#'                       ras = rasterL,
#'                       corTEST = 0.7,
#'                       vifTEST = TRUE,
#'                       corVIF = 0.7,
#'                       use = "complete.obs")
#'
#' @export
wsl.varKeep=function(PPower.object,ras,corTEST=0.7,vifTEST=FALSE,corVIF=0.7,...)
{
  # Looping over the PPower.object output of the varPPower function

  perVAR=list()

  for (i in 1:length(PPower.object))
  {
    cat ("###### Class object","'",names(PPower.object)[i],"'","\n")

    # Focus on the summary line of our object and ordering

    Q1a=row.names(PPower.object[[i]]) %in% "meanSummary"
    t1a=as.numeric(PPower.object[[i]][Q1a,])

    Q1b=row.names(PPower.object[[i]]) %in% "sdSummary"
    t1b=as.numeric(PPower.object[[i]][Q1b,])

    t2=data.frame(PP=t1a,SD=t1b,VAR=colnames(PPower.object[[i]]),stringsAsFactors=F)
    t3=t2[order(t2$PP,decreasing=T),]

    # Create duplicate

    tref=t3
    tref$rank=1:nrow(tref)

    # Extract raster name(s)

    if (class(ras) %in% "list")
    {
      rasID=sapply(ras,function(x) names(x))
    } else {
      rasID=names(ras)
    }

    # Test descending colinearity between ranked predictors

    toKeep=list()
    toKeepNot=list()

    for (j in 1:length(t2$VAR))
    {
      cat ("cor.test 'selection-deletion' n°",j,"out of",length(t2$VAR),"\n")

      varA=ras[[which(rasID %in% t3[1,"VAR"])]]
      varB=ras[[which(rasID %in% t3[2,"VAR"])]]
      # Maybe a faster cor() needed
      test=cor(values(varA),values(varB),use="complete.obs")

      if (j!=(length(t2$VAR)-1))
      {
        # Argument for selecting-deleting predictors
        if (abs(test)>corTEST) {
          toKeepNot[[j]]=names(varB)
          t3=t3[!(t3$VAR %in% names(varB)),]
        } else {
          toKeep[[j]]=names(varA)
          t3=t3[!(t3$VAR %in% names(varA)),]
        }
      } else {
        # Argument for last comparison
        if (abs(test)>corTEST) {
          toKeepNot[[j]]=names(varA)
        } else {
          toKeep[[j]]=c(names(varA),names(varB))
        }
      }
    }
    # Create a vectors from lists
    toKeepNot=unique(unlist(toKeepNot))
    toKeep=unique(unlist(toKeep))

    # Keep only the non-NAs
    toKeepNot=suppressWarnings(toKeepNot[!is.na(toKeepNot)])
    toKeep=suppressWarnings(toKeep[!is.na(toKeep)])

    # Generates data.frames
    r1=data.frame(VAR=toKeep,stringsAsFactors=F)
    r1=try(merge(r1,tref,by="VAR",sort=FALSE,all.x=T),silent=T)
    if (class(r1) %in% "try-error"){r1=NULL}

    r2=data.frame(VAR=toKeepNot,stringsAsFactors=F)
    r2=try(merge(r2,tref,by="VAR",sort=FALSE,all.x=T),silent=T)
    if (class(r2) %in% "try-error"){r2=NULL}

    # Add description titles
    toBeOrNot=list(r1,r2)
    names(toBeOrNot)=c("no_correlation_PPrank","correlation_PPrank")

    # Store
    if (!vifTEST) {
      perVAR[[i]]=toBeOrNot
    } else {

      cat("vifTEST=TRUE","\n")

      K=toBeOrNot$no_correlation_PPrank$VAR
      rasK=sapply(ras,function(x) names(x) %in% K)
      stacK=stack(ras[rasK])
      toBeOrNot$VIFselect=vifcor(stacK,corVIF)
      perVAR[[i]]=toBeOrNot
    }
  }

  # Add predictor classes description
  names(perVAR)=names(PPower.object)

  return(perVAR)
}



