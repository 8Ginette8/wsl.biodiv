# Libraries

library(raster)
library(sp)
library(ecospat)
library(foreach)
library(doParallel)
library(usdm)

###############################
##### First function
###############################

# Description: 'varPPower' is to obtain a ranking of predictive power
# (R2 or D2) of predictors (or class) trying to explain species
# observations (or other ecological data such diversity-count or
# continuous variables)

# "object" = vector() for matrix/df of xy coords
# "object" = list() for SpatialPoints
# points, species, val, weight should be of same length

# points    	=====> matrix/df of xy coords; SpatialPoints; SpatialPoints list
# val       	=====> "object" of presences-absences, or of other spatial infos
# species   	=====> vector of species names or other infos
# ras 		  =====> one RasterLayer or a list of RasterLayer
# rasCLASS	=====> a vector associating a Class ID to the raster(s)
# mlinear	  =====> using lm() or glm() depending on the type of data
# glmMODE	  =====> choosing the glm link function: "binomial" (defaut) or "poisson"
# weight  	=====> "object" of weight information associated with "val" to use in the model fit
# poly 	 	  =====> using the polynomial quadratic function if we want a more flexible fit
# polyBRUT  	=====> Force poly function when NAs are found
# cores     	=====> number of cores for parallelisation of extract()
# ...		    =====> additional arguments to add to lm() or glm() function

### ==================================================================
### Loop over each species to extract the D2/R2 of each predictor
### ==================================================================

wsl.varPPower=function(points,val,species=NULL,ras,rasCLASS=NULL,
	mlinear=FALSE,glmMODE="binomial",weight=NULL,poly=FALSE,
	polyBRUT=FALSE,parallel=FALSE,cores=detectCores()/2,...)
{
	### ==================================================================
	### Error handlings
	### ==================================================================

	# Checking class of input rasters
	# Checking if raster have the same CRS
	if (!(class(ras) %in% "list")) {
		cond1=class(ras)[1] %in% "RasterLayer"
		cond2=TRUE
	} else {
		r.class=unlist(lapply(ras,function(x) class(x)))
		cond1=all(r.class=="RasterLayer")
		r.crs=unlist(lapply(ras,function(x) proj4string(x)))
		cond2=length(unique(r.crs))==1
	}

	if (cond1==FALSE) {
		stop("Input raster(s) should be of class 'RasterLayer' !")
	}
	if (cond2==FALSE) {
		stop("Input raster(s) should have the same CRS !")
	}

	# Checking if input xy are correct
	cond3=!(class(points)[1] %in% "SpatialPoints")
	cond4=!(ncol(points)==2 || all(colnames(points)%in%c("x","y")))

	if (cond3 & cond4) {
    	stop("Supplied points should be of class 'SpatialPoints' or a data.frame/matrix with two columns named x and y!")
    }

    # Checking if raster and points have the same CRS
    if (any(class(points) %in% c("matrix","data.frame"))) {
    	if (!(compareCRS(ras[[1]],CRS("+init=epsg:4326")))) {
    		stop("Input raster(s) and XY points should be in Longitude/Latitude.\nPlease use the following CRS:'+init=epsg:4326'")
    	}
    } else {
    	cond5=c(class(points) %in% "list",class(ras) %in% "list")

    	crsP=ifelse(cond5[1],unique(sapply(points,function(x) proj4string(x))),proj4string(points))
    	crsR=ifelse(cond5[2],unique(sapply(ras,function(x) proj4string(x))),proj4string(ras))
    	
    	if(!(compareCRS(crsP,crsR))) {
			stop("Input raster(s) and supplied 'SpatialPoints' should have the same CRS")
		}
    }

	### ==================================================================
	### Conditions for 'species' and 'rasCLASS' parameters
	### ==================================================================

	# Checking if the "species" parameter is not NULL
	if (is.null(species)) {
		spN=1
	} else {
		
		# Checking Length of 'points' and 'species' parameters
		spN=unique(species)
		if (length(spN)==1) {
			stop("'species' parameter should be NULL or have more than one species...")
		}

		LL=sapply(list(species,val,weight),function(x) length(x))
		
		if (any(class(points) %in% c("matrix","data.frame"))) {
			if (!(all(LL==nrow(points)))) {
				stop("XY points, 'species', 'val' & 'weight' should be of same length...")
			}
		} else {
			if (!(all(LL==length(points)))) {
				stop("SpatialPoints, 'species', 'val' & 'weight' should be of same length...")
			}
		}
	}

	# Checking if the "r.class" parameter is not NULL and assigning classes
	if (is.null(rasCLASS)) {
		VARclass="No class"
	} else {
		
		# Checking length of 'ras' and 'r.class' parameters
		VARclass=unique(rasCLASS)
		if (length(VARclass)==1) {
			stop("'rasCLASS' parameter should NULL or have more than one class...")
		}
		if (!(length(rasCLASS)==length(ras))) {
			stop("'rasCLASS' & 'ras' parameters should be of same length...")
		}
	}

	### ==================================================================
	### Infos parallelisation
	### ==================================================================

	if (parallel)
	{
		cat("Parallelisation is TRUE with",cores,"cores","\n")
	} else {
		cat("Parallelisation is FALSE","\n")
	}

	### ==================================================================
	### Applying models
	### ==================================================================

	D1VAR=list()

	# Looping over the raster classes
	for (i in 1:length(VARclass))
	{
		cat("####### Processing class ","'",VARclass[i],"'",
			" (",i," out of ",length(VARclass)," Classes)...",sep="","\n")

	# Extract rasters from class i
		if (is.null(rasCLASS)) {
			VARn=ras
		} else {
			Q2=which(rasCLASS %in% VARclass[i])
			VARn=ras[Q2]
		}

		D2SP=list()

		# Looping over the species
		for (j in 1:length(spN))
		{
			cat("Processing species ","'",unique(species)[j],"'",
				" (",j," out of ",length(spN)," Species)...",sep="","\n")

			# Keeping for each species the right info
				# If no species infos
			if (is.null(species)) {
				param=list(points,val,weight)
			} else {
				Q1=which(species %in% spN[j])

				# If points are XY coords
				if (any(class(points) %in% c("matrix","data.frame"))) {
					param=list(points[Q1,],val[Q1],weight[Q1])

				# If points are SpatialPoints
				} else {
					param=list(points[Q1][[1]],val[Q1][[1]],weight[Q1][[1]])
				}
			}

			# Extracting for every XY the values of all rasters from the class

				# Choosing to go Parallel
				if (class(VARn) %in% "list") {
					# TRUE
					if (parallel)
					{
						# Creating a cluster with cpus
						cl <- makeCluster(cores)
						registerDoParallel(cl)

						# Going Parallel
						val.X=foreach (k=1:length(VARn),.packages=c("doParallel","foreach","raster")) %dopar%
						{
							out=extract(VARn[[k]],param[[1]],method="simple")
							return(out)
						}

						# Closing the cluster
						stopCluster(cl)

					# FALSE	
					} else {
						val.X=foreach (k=1:length(VARn)) %do%
						{
							cat("extract",names(VARn[[k]]),"\n")
							
							out=extract(VARn[[k]],param[[1]],method="simple")
							return(out)
						}
					}

				# Simple extract in case of just one raster
				} else {
					val.X=list(extract(VARn[[1]],param[[1]],method="simple"))
				}

				# Assign names
			names(val.X)=sprintf("X%d",1:length(val.X))

			# Write formula for each raster val.X
			inputName=paste0("val.X$",names(val.X))

			if (poly) {
				if (!polyBRUT)
				{
					form.m=lapply(inputName,function(x) as.formula(paste("param[[2]]",paste("poly(",x,",2)"),sep="~")))
				} else {
					form.m=lapply(inputName,function(x) as.formula(paste("param[[2]]",paste(x,"+ I(",x,"^2)"),sep="~")))
				}
			} else {
				form.m=lapply(inputName,function(x) as.formula(paste("param[[2]]",paste(x),sep="~")))
			}
			
			# Apply formula GLM or LM
			if (mlinear) {
				cat("Calculating R2 adj.lm...",sep="","\n")
				DR2=sapply(form.m,function(x) summary(lm(x,weights=param[[3]],...))$adj.r.squared)
			} else {
				cat("Calculating D2 adj.glm...",sep="","\n")
				DR2=sapply(form.m,function(x) ecospat.adj.D2.glm(glm(x,weights=param[[3]],family=glmMODE,...)))*100
			}

			# Store
			D2SP[[j]]=DR2
		}

		# Aggregating list and assigning infos
		names(D2SP)=unique(species)
		result=do.call("rbind",D2SP)

		# Rename columns
		if (!(ncol(result)==1)) {
			colnames(result)=sapply(VARn,function(x) names(x))
		} else {
			colnames(result)=names(VARn)
		}
		
		# Do the summary
		if (nrow(result)==1) {
			meanSummary=NULL
			sdSummary=NULL
			row.names(result)="meanSummary"
		} else {
			# Creating a new row with the PPower summary
			meanSummary=apply(result,2,mean)
			sdSummary=apply(result,2,sd)
		}

		# Store species result for each assigned class of predictors
		result2=rbind(result,meanSummary,sdSummary)
		D1VAR[[i]]=result2
	}
	names(D1VAR)=VARclass

	# Return the results
	return(D1VAR)
}


###############################
##### Second function
###############################

# Description: 'varKeep' is a function to select the predictors
# used to explain the data distribution in 'varPPower', according
# to their mean PP ranking and their successive descending correlation.

# var.object =====> output of the varPPower function
# ras        =====> same raster or list of rasters used in 'varPPower'
# corTEST	   =====> level of accepted correlation to select our predictors
# ...		     =====> additional arguments to add to the cor() function

wsl.varKeep=function(PPower.object,ras,corTEST=0.7,vifTEST=NULL,...)
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

		for (j in 1:(length(t2$VAR)-1))
		{
			cat ("cor.test 'selection-deletion' n°",j,"out of",(length(t2$VAR)-1),"\n")

			varA=ras[[which(rasID %in% t3[1,"VAR"])]]
			varB=ras[[which(rasID %in% t3[2,"VAR"])]]
			# Maybe a faster cor() needed
			test=cor(values(varA),values(varB),...)

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
		toKeepNot=unlist(toKeepNot)
		toKeep=unlist(toKeep)

		# Keep only the non-NAs
		toKeepNot=toKeepNot[!is.na(toKeepNot)]
		toKeep=toKeep[!is.na(toKeep)]

		# Generates data.frames
		r1=data.frame(VAR=toKeep,stringsAsFactors=F)
		r1=merge(r1,tref,by="VAR",sort=FALSE,all.x=T)

		r2=data.frame(VAR=toKeepNot,stringsAsFactors=F)
		r2=merge(r2,tref,by="VAR",sort=FALSE,all.x=T)
		
		# Add description titles
		toBeOrNot=list(r1,r2)
		names(toBeOrNot)=c("no_correlation_PPrank","correlation_PPrank")

		# Store
		if (!vifTEST) {
			perVAR[[i]]=toBeOrNot
		} else {

			cat("vifTEST=TRUE","\n")

			K=toBeOrNot$no_correlation_PPrank$VAR
			rasK=sapply(rasterL,function(x) names(x) %in% K)
			stacK=stack(rasterL[rasK])
			toBeOrNot$VIFselect=vifstep(stacK)
			perVAR[[i]]=toBeOrNot
		}
	}

	# Add predictor classes description
	names(perVAR)=names(PPower.object)
	
	return(perVAR)
}
