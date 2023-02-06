### =========================================================================
### wsl.obs.filter
### =========================================================================
#' Filter set of observations
#'
#' Filter presences / absences through a chosen resolution grid.
#'
#' @param o.xy Object of class 'matrix' or 'data frame' with two columns named
#' "x" and "y". Must be used alone to apply grid filetring to one set of observations
#' (presences or absences).
#' @param a.xy Object of class 'matrix' or 'data frame' with two columns named "x" and
#' "y". NULL by default. If used, "o.xy" becomes presences and "a.xy" absences.
#' Resolution grid filtering is applied to both observation sets, but coordinates of
#' "a.xy" falling within the same grid cells as "o.xy" are removed.
#' @param grid Object of class 'RasterLayer', 'RasterBrick' or 'RasterStack' of
#' desired resolution and extent.
#' @return Object of class 'matrix' or 'data frame' with two columns named "x" and "y"
#' comprising the new set of observations filtered at grid resolution. If "a.xy" is
#' used the output is a list of two data.frame: filtered presences and filtered absences
#' @author Yohann Chauvier
#' @examples
#' 
#' ### Load my binary observations species data
#' 
#' library(raster)
#' 
#' data(var_select_XYtest)
#' data(exrst)
#' 
#' ### wsl.obs.filter(): example for the first species
#'
#'    # Loading observations: presences and absences
#' 
#' presences = coordinates(mySP[[1]])[myPA[[1]] %in% "1",]
#' absences = coordinates(mySP[[1]])[myPA[[1]] %in% "0",]
#'
#'    # Loading grid
#' 
#' r.layer = rst[[1]]
#'
#'    # To filter observations by the grid only
#' 
#' pres.filtered = wsl.obs.filter(presences,grid=r.layer)
#' abs.filtered = wsl.obs.filter(absences,grid=r.layer)
#'
#'    # To filter observations by the grid & remove abs. in cells where we also find pres.
#' 
#' PresAbs.filtered = wsl.obs.filter(presences,absences,r.layer)
#'
#'    # Count presences (same filtering)
#' 
#' nrow(PresAbs.filtered[[1]])
#' nrow(pres.filtered)
#'
#'    # Count Absences (filtering plus removal of duplicated absences)
#' 
#' nrow(abs.filtered)
#' nrow(PresAbs.filtered[[2]])
#'
#'    # Visual
#' 
#' par(mfrow=c(1,2))
#' plot(presences)
#' plot(pres.filtered)
#'
#' par(mfrow=c(1,2))
#' plot(absences)
#' plot(abs.filtered)
#' 
#' @export
wsl.obs.filter=function(o.xy,a.xy=NULL,grid)
{
  # Check 'o.xy' input

  if(ncol(o.xy)!=2 || !all(colnames(o.xy)%in%c("x","y"))){
    stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
  }

    # Check 'ras' input

    if(!(class(grid)%in%c("RasterBrick","RasterStack","RasterLayer"))) {
      stop("env.layer should be of class RasterBrick, RasterStack or RasterLayer!")
    }

    # Apply simple filtering
    if (!is.null(a.xy)) {
      
      # Check 'a.xy' input

      if(ncol(a.xy)!=2 || !all(colnames(a.xy)%in%c("x","y"))){
        stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
      }

      # Position presences and absences
      posP=cellFromXY(grid,o.xy)
      posA=cellFromXY(grid,a.xy)

      # Remove absences where we find presences
      posA=posA[!(posA %in% posP)]

      # Extract new presences/absences and regroup
      new.pxy=coordinates(grid)[unique(posP),]
      new.axy=coordinates(grid)[unique(posA),]
      new.oxy=list(new.pxy,new.axy)
      names(new.oxy)=c("Presences","Absences")

    } else {
      # Apply simple filtering
      posCELL=cellFromXY(grid,o.xy)
      new.oxy=coordinates(grid)[unique(posCELL),]
    }

    if (class(new.oxy)[1]%in%"numeric"){
      new.oxy=matrix(new.oxy,ncol=2)
      colnames(new.oxy)=c("x","y")
    }

    return(new.oxy)
}
