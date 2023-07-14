### =========================================================================
### get_env
### =========================================================================
#' Extract environments based on yearly observations (local or URL sources)
#'
#' Based on the integrated list of CHELSA files (https://chelsa-climate.org/;
#' CHELSAcruts), this function extracts environmental values from given observations
#' at their year of sampling (mean, median, min, max or montly averaged).It is
#' important to note that CHELSAcruts is available for 1901-2016.
#'
#' @param env_layers Vector. Which environment should be extracted? One or several may be
#' chosen among 'prec', 'tmax' and 'tmin'. Default is c('prec','tmin').
#' @param timeline String. What type of data extraction is applied? Available choices are
#' 'yearly' or 'monthly'. One or the other will give mean, median, min and max values. To
#' note that 'monthly' will run faster than 'yearly' since less files will be downloaded.
#' Default is 'monthly'.
#' @param months Vector. Only if time = 'monthly'. Which months should be accounted
#' for? Input must be from 1 to 12. To note that values from -1200 (December n-1) to 1200
#' (December n+1; step of 100) may be chosen, and corresponds to the months of the year n-1
#' and n+1 respectively. Default is June, July and August.
#' @param obs A data.frame with three columns (describing in order: Longitude, Latitude and year
#' of sampling). The observations needs to be from the same CRS as the 'env_layers'.
#' @param download Logical. If FALSE, the environmental layers are locally available and will be
#' extracted from 'd.path'.
#' @param d.path Character. The folder path where the environmental layers should be saved or stored.
#' @return A data.frame describing for each observations their associated environmental values.
#' @references 
#' ...
#' @examples
#' 
#' # Packages
#' library(wsl.biodiv)
#' data(xy_ppm)
#'  
#' # Dummy observations
#' random.year = (1901:1921)[sample(1:21,nrow(xy.ppm),replace=TRUE)]
#' obs = data.frame(xy.ppm[,c("x","y")],year=random.year)[1:5,]
#' 
#' # Run function
#' sp.env = get_env(env_layers = 'tmin',
#'                  timeline = 'monthly',
#'                  months = c(-1200,1,2),
#'                  obs = obs,
#' 					download = TRUE,
#' 					d.path = getwd())
#' 
#' @export
#' 
get_env = function(env_layers = c('prec',"tmin"),
				   timeline = 'monthly',
				   months = c(6,7,8),
				   obs = data.frame(),
				   download = TRUE,
				   d.path = getwd())
{
	# Inventory of CHELAcruts files
	url.base = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/"
	cf = expand.grid("CHELSAcruts",c("prec","tmin","tmax"),1:12,1901:2016,"V.1.0.tif")
	cruts.files = paste(cf$Var1,cf$Var2,cf$Var3,cf$Var4,cf$Var5,sep="_")
	cruts.files = paste0(c("prec/","tmin/","tmax/"),cruts.files)
	cruts.files = paste0(url.base,cruts.files)

	# In case the year before or after is requested
	base.year = unique(obs$year)
	if (any(months >= 100)&any(months < 0)) {
		stop("Average with months n+1 & n-1 not allowed...")#add.year = unique(c(obs$year-1,obs$year+1))
	} else if (any(months >= 100)) {
		add.year = unique(obs$year)+1
	} else if (any(months < 0)) {
		add.year = unique(obs$year)-1
	}

	# Extract the files we want
	tar.env = cruts.files[grepl(paste(env_layers,collapse="|"),cruts.files)]
	tar.year = unlist(lapply(base.year,function(x) tar.env[grepl(x,tar.env)]))
	if (timeline == "monthly") {
		tar.month = unlist(lapply(months,function(x) tar.year[grepl(paste0("_",x,"_"),tar.year)]))
		if (exists("add.year")) {
			tar.add = unlist(lapply(add.year,function(x) tar.env[grepl(x,tar.env)]))
			months.n = abs(months[months >= 100|months < 0]/100)
			tar.mn = unlist(lapply(months.n,function(x) tar.add[grepl(paste0("_",x,"_"),tar.add)]))
			tar.year = unique(c(tar.month,tar.mn))
		} else {
			tar.year = unique(tar.month)
		}
	}

	# URL or files
	if (download) {

		# Main path to download and/or extract the rasters
		out.path = paste0(d.path,"/",gsub(".*/","",tar.year))

		# Loop
		lapply(1:length(tar.year), function(x){

			# Download
			ftp.file = try(download.file(tar.year[x],out.path[x],mode="wb"))
			if (is(ftp.file,"try-error"))
			{
				while(is(ftp.file,"try-error"))
				{
					cat("Reconnecting...","\n",sep="")
					ftp.file = try(download.file(tar.year[x],out.path,mode="wb"))
				}
			}
		})

	} else {
		out.path = list.files(d.path)
	}

	# Filter the files according to the parameters
	out.path = out.path[grepl(".tif",out.path)]
	if (exists("months.n")) {e.k = unique(c(months,months.n))} else {e.k = unique(months)}
	out.path = out.path[unlist(sapply(e.k,function(x) grep(paste0("_",x,"_"),out.path)))]

	# Open the rasters
	our.ras = lapply(out.path, function(x) rast(x))

	# Prepare out data.frame
	out.df = as.data.frame(matrix(nrow=nrow(obs),ncol=ncol(obs)+length(env_layers)*4))
	out.df[,1:3] = obs
	egg = expand.grid(c("mean","median","min","max"),env_layers)
	names(out.df) = c(names(obs),sprintf('%s.%s',egg[,2],egg[,1]))

	# Years that needs extraction
	ext.years = unique(out.df[,3])

	# Loop over environmental variable
	pos.col = c(0,3,6)
	for (i in 1:length(env_layers)) {

		# Get the right rasters
		list.p1 = gsub(".*/","",sapply(our.ras,sources))
		id.ras = our.ras[grepl(env_layers[i],list.p1)]
		list.p2 = list.p1[grepl(env_layers[i],list.p1)]

		# Loop over years
		for (j in 1:length(ext.years)) {

			# Observations and predictors from this year (or previous/next)
			year.obs = out.df[out.df[,3]%in%ext.years[j],1:2]
			year.ras = try(rast(id.ras[grepl(ext.years[j],list.p2)]),silent=TRUE)
			if (exists("add.year")) {

				# Extract environment from previous/next year
				if (any(months<0)) {
					months.ras = rast(id.ras[grepl(ext.years[j]-1,list.p2)])
				} else {
					months.ras = rast(id.ras[grepl(ext.years[j]+1,list.p2)])
				}

				# Keep 
				toK = sapply(months.n,function(x) grep(paste0("_",x,"_"),names(months.ras)))
				months.rasK = months.ras[[toK]]

				# Remove not corresponding months or not
				if (class(year.ras)%in%"try-error"){
					year.ras = months.rasK
				} else {
					toR = sapply(months.n,function(x) grep(paste0("_",x,"_"),names(year.ras)))
					q1 = abs(months[months>=100|months<0]/100)
					q2 = months[!months>=100|months<0]
					if (!class(toR)%in%"list") {year.ras = c(year.ras,months.rasK)[[-toR]]}
					if (class(toR)%in%"list"|q1%in%q2) {year.ras = c(year.ras,months.rasK)}
				}
			}

			# Extract infos
			env.info = extract(year.ras,year.obs)
			env.val = as.numeric(env.info[,!names(env.info)%in%"ID"])
			out.df[out.df[,3]%in%ext.years[j],3+i+pos.col[i]] = mean(env.val,na.rm=TRUE)
			out.df[out.df[,3]%in%ext.years[j],3+i+1+pos.col[i]] = median(env.val,na.rm=TRUE)
			out.df[out.df[,3]%in%ext.years[j],3+i+2+pos.col[i]] = min(env.val,na.rm=TRUE)
			out.df[out.df[,3]%in%ext.years[j],3+i+3+pos.col[i]] = max(env.val,na.rm=TRUE)
		}
	}

	# Return
	return(out.df)
}

