### =========================================================================
### Environmental Bias Correction (EBC)
### =========================================================================
#' Correct environmental bias of species observations via environmental clustering
#'
#' Our wsl.ebc function corrects environmental bias from an observational dataset based on
#' environmental stratification / clustering. A map of n clusters is generated based on input
#' raster layers (in the article, climate predictors used in PPMs are employed as inputs).
#' Following random equal-stratified sampling design (EBCe; Hirzel & Guisan, 2002), for each
#' cluster, the number of observations per species (relative to all others) is artificially
#' rescaled to the total number of observations found in the cluster presenting the highest
#' observation density. Proportional-stratified sampling design (EBCp; Hirzel & Guisan, 2002)
#' adds a second step: for each cluster, the number of observations per species may additionally
#' be multiplied by the cluster's area (i.e. proportion of pixels in percentage relative to the
#' study area), or by its logarithm (default; consensus between EBCe and EBCp). Resulting output
#' indicates a new number of observations per cluster and species, that the function automatically
#' sub-samples with replacement over the original observational dataset. This function may be used
#' for presences and absences distinctively.
#'
#' @param obs Object of class matrix or data frame with three columns named sp.id (character),
#' x (numeric) and y (numeric). More than one observation per species must be referenced.
#' @param ras Object of class RasterBrick, RasterStack, list of RasterLayer of desired resolution
#' and extent. Used to generate the map of clusters needed to summarize the environmental space of
#' the study area.
#' @param pportional Logical. Should environmental stratification of observations be proportional
#' to the clusters' areas? If TRUE, EBCp applies.
#' @param plog Logical. Should EBCp apply with a logarithm? If TRUE, a stratification consensus
#' between EBCp and EBCe applies.
#' @param nclust Number of chosen clusters. Default is 50.
#' @param sp.specific Logical. Should EBC apply only for species whose environmental bias follows
#' the overall one (i.e. the number of original species observations per cluster is correlated
#' with that of the full dataset)?
#' @param sp.cor  If sp.specific = TRUE, spearman's correlation tests are by default set to 0.5;
#' i.e. species with r < 0.5 are excluded from the function outputs.
#' @param keep.bias Default is TRUE. Strongly recommended to use when sp.cor = TRUE. Per species,
#' should the number of observations of the cluster in which the species was originally sampled
#' most often, be preserved? Said differently, for each species, should the cluster with the most
#' original species observations be as representative as the cluster with the most corrected
#' observations? If TRUE, after EBC applies, the number of observations per species in their
#' densest original cluster is set to that of the densest corrected cluster.
#' @param filter Logical. Should the observations be filtered according to 'ras' resolution?
#' @param path Path folder where the new species observation files should be saved.
#' @param ... Additional arguments passed on to the ‘clara’ function (package 'cluster') 
#' @return The function returns in ‘path’, one text file of corrected observations (presences
#' or absences) per species. If the number of new EBC observations per species is too large,
#' sampling those randomly without replacement before model calibrations is advised.
#' @author Yohann Chauvier
#' @examples
#' 
# # Cleaning
# rm(list = ls()); graphics.off()

# # max.print
# options(max.print=500)

# # Load libraries
# library(cluster)
# library(raster)

# # Load data
# load("exrst.RData")
# load("xy_ppm.RData")
# names(xy.ppm)[1] = "sp.id"

# # Run EBC function with the log consensus
# # sp.specific = TRUE --> Select corrected XY outputs in a species-specific manner
# # keep.bias = TRUE   --> Preserve initial observer bias of each species
# wsl.ebc(obs = xy.ppm,
#  		ras = rst[[1:5]],
#  		pportional = TRUE,
#  		plog = TRUE,
#  		nclust = 50,
#  		sp.specific = TRUE,
#  		sp.cor = 0.5,
#  		filter = TRUE,
#  		keep.bias = TRUE,
#  		path = getwd())
 
# # Open corrected observations
# files = list.files(getwd())
# target.files = files[grep("_obs_corrected_",files)]
# correct.obs = lapply(target.files, function(x) obs=read.table(paste0(getwd(),"/",x)))
# correct.obs = do.call("rbind",correct.obs)

# # Extract with corrected and uncorrected observations environmental values of 'rst'
# correct.env = extract(rst,correct.obs[,-1])
# uncorrect.env = extract(rst,xy.ppm[,-1])

# # Sampling of the environment (e.g. second layer) and of environmental values extracted with observations
# f.correct.env = correct.env[,2]
# f.uncorrect.env = uncorrect.env[,2]
# first.env = rst[[1]][!is.na(rst[[2]])]
# samp.env = first.env[sample(1:length(first.env),2e3)]
# samp.vc=f.correct.env[sample(1:length(f.correct.env),2e3)]
# samp.vnc=f.uncorrect.env[sample(1:length(f.uncorrect.env),2e3)]
 
# # Plot distribution histograms for each 
# par(mfrow=c(2,2))
# hist(samp.env,breaks=10)
# hist(samp.vc,breaks=10)
# hist(samp.vnc,breaks=10)
#' 
#' @export
wsl.ebc=function(obs=NULL,ras=NULL,pportional=TRUE,plog=TRUE,nclust=50,
	sp.specific=TRUE,sp.cor=0.5,keep.bias=TRUE,filter=FALSE,path=NULL,...)
{
	# 'nclust' must be at least equal to '2'
	if (nclust<=1 || is.null(nclust)){
		stop("'nclust' must be equal to 2 or more...!")
	}

	# 'obs' must be a data.frame of observations with column "sp.id","x" and "y"
	if (!(class(obs)%in%c("data.frame","matrix"))){
		stop("'obs' must be an object of class 'data.frame' or 'matrix'...!")
	}

	if (class(obs)%in%"matrix") {
		obs=as.data.frame(obs,stringsAsFactors=FALSE)
		obs$x=as.numeric(obs$x)
		obs$y=as.numeric(obs$y)
	}
	
	if (any(table(obs$sp.id)%in%1)){
		stop("Each species must have more than one observation...!")
	}

	# 'ras' must be of 'brick' or 'stack' class or list of 'rasters'
	if (!(class(ras)%in%c("RasterBrick","RasterStack","list"))){
		stop("'ras' must be an object of class 'RasterBrick', 'RasterStack' or 'list'...!")
	} else if (length(ras)%in%1 || nlayers(ras)%in%1){
		stop("'ras' must include more than one 'RasterLayer'...!")
	}
	if (class(ras)%in%"list"){
		ras=stack(ras)
	}

	# First info message
	cat("Step 1 --> CLARA algorithm processing...","\n")
	
	# Convert raster in the right CLARA format
	id.toReplace=complete.cases(ras[])
	toClara=ras[][id.toReplace,]

	# Run CLARA and assign results to right pixels
	blocks=clara(toClara,nclust,...)
	toNew.ras=ras[[1]]
	toNew.ras[][id.toReplace]=blocks$clustering
	toNew.ras[][!id.toReplace]=NA

	# Start inventorying species observations per cluster 
	sp.ref=as.character(unique(obs$sp.id))

	# Second info message
	cat("Step 2 --> Species inventory...","\n")

	cluster.infos=list()
	cluster.tab=list()
	for (i in 1:length(sp.ref))
	{
		# Extract each species one by one & filter according to grid resolution if TRUE
		sp.target=obs[obs$sp.id%in%sp.ref[i],c("x","y")]
		if (filter){
			sp.target=wsl.obs.filter(sp.target[,c("x","y")],grid=toNew.ras)
		}

		# Extract cluster values with warnings if XY outside raster extent or assign to NA values
		clustV=extract(toNew.ras,sp.target)
		if(any(is.na(clustV))) {
			warning(paste("XY coordinates of",sp.ref[i],
				"outside 'ras' boundaries, or have extracted NAs environmental values...",sep=" "))
		}
		summaryV=table(clustV)

		# Fill a vector of length 'nclust'
		sp.vec=rep(0,nclust)
		sp.vec[1:nclust%in%as.numeric(names(summaryV))]=summaryV

		# Store infos
		cluster.infos[[i]]=cbind(sp.target,clustV)
		cluster.tab[[i]]=c(sp.ref[i],sp.vec)
	}

	# Combine and apply a sum to all columns
	cluster.all=as.data.frame(do.call("rbind",cluster.tab))
	cluster.all[,2:(nclust+1)]=sapply(cluster.all[,2:(nclust+1)],
		function(x) as.numeric(as.character(x)))
	cluster.sum=apply(cluster.all[,-1],2,sum)

	# Which cluster is max?
	plateau=max(cluster.sum)

	# Calculate area values based on n pixels per cluster (formated)
	pixels.clust=table(toNew.ras[])
	area.clust=pixels.clust/sum(pixels.clust)*100

	# Get reduction coef for each cluster specific to "plateau"
	reducoef=cluster.sum/plateau

	# Apply it to the species table and round up to the integer
	down.l=lapply(1:length(reducoef),function(x) cluster.all[,1+x]/reducoef[x])
	cluster.down=do.call("cbind",down.l)
	cluster.down=ceiling(cluster.down)

	# Weight observations by cluster areas
	if (pportional){
		cat ("...Environmental proportional stratification = TRUE...","\n")

		if (plog){
			cat ("...log...","\n")
			area.clust=log(area.clust+1)
		}

		down.l2=lapply(1:length(reducoef),function(x) cluster.down[,x]*area.clust[x])
		cluster.down2=do.call("cbind",down.l2)
		cluster.down=ceiling(cluster.down2)

	} else {
		cat ("...Environmental proportional stratification = FALSE...","\n")
	}

	# Replace NAs by 0 in case absolutely no informations fall in one cluster
	cluster.down[is.na(cluster.down)]=0

	# Store
	cl.out=data.frame(x0=cluster.all[,1],cluster.down)
	cl.out[,1]=as.character(cl.out[,1])

	# Do we apply EBC only for target species?
	if (!(sp.specific)){
		cat("Step 3 --> EBC set to non species-specific...","\n")
	} else {
		cat("Step 3 --> EBC set to species-specific...","\n")

		# Tchao species whose environmental bias is not correlated with the general one
		spcl.temp=
		lapply(1:nrow(cluster.all),function(x) {
			sp.sumo=as.numeric(cluster.all[x,-1])
			cor.ebc=cor(sp.sumo,cluster.sum)
			if (cor.ebc>=sp.cor | is.na(cor.ebc)){
				return(cl.out[x,])
			}
		})
		# Format cluster.infos & bind results
		ebc.null=sapply(spcl.temp,is.null)
		cluster.infos=cluster.infos[!ebc.null]
		cl.out=do.call("rbind",spcl.temp)
	}

	# Should the initial environmental bias be preserved?
	if (!keep.bias){
		cat("Step 4 --> Initial environmental bias per species is not preserved...","\n")
	} else {
		cat("Step 4 --> Initial environmental bias per species is preserved...","\n")

		# Keeping initial biased environment
		spcl.temp2=
		lapply(1:nrow(cl.out),function(x) {
			spcl.out=cl.out[x,]
			o.clust=cluster.all[cluster.all[,1]%in%spcl.out[,1],]
			maxo=which.max(o.clust[,-1])[1]
			spcl.out[,maxo+1]=max(spcl.out[,-1],na.rm=TRUE)
			return(spcl.out)
		})
		# Format cluster.infos & bind results
		cl.out=do.call("rbind",spcl.temp2)
	}

	# Fourth infos message
	cat("Step 5 --> Environmental Bias Correction (EBC) starting...","\n")
	
	# Write new observation files per species in choosen 'path'
	for (i in 1:nrow(cl.out))
	{
		# Extract number of obs. per cluster to sample
		sp.infos=cl.out[i,]
		sp.save=as.character(sp.infos[1])

		# Extract original observations from species
		cl.i=cluster.infos[[i]]

		# Fourth infos message
		cat("Processing ",sp.save,"...","\n")

		# Extract obs. according to each cluster
		Cselect=list()
		for (j in 1:nclust) {

			# Extract observations in clusters
			tg=cl.i[cl.i[,"clustV"]%in%j,c("x","y"),drop=FALSE]

			# Sample with replacements
			n.samp=as.numeric(sp.infos[-1][j])
			do.samp=tg[sample(1:nrow(tg),n.samp,replace=TRUE),c("x","y")]

			# Store with cluster ID depending on the output
			if (nrow(do.samp)%in%0){
				next
			} else if (class(do.samp) %in% "numeric") {
        		do.samp=matrix(do.samp,ncol=2)
        		colnames(do.samp)=c("x", "y")
    		}

			Cselect[[j]]=data.frame(cluster.id=paste0("clust.",j),do.samp,stringsAsFactors=FALSE)
		}

		# Rbind & save
		new.obs=do.call("rbind",Cselect)
		out.sp=paste0(path,"/",Sys.Date(),"_obs_corrected_",sp.save,".txt")
		write.table(new.obs,out.sp)
	}
}