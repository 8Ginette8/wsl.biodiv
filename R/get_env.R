### =========================================================================
### get_env
### =========================================================================
#' Extract environments based on yearly observations (local or URL sources)
#'
#' Based on local files or URL sources, this function extracts environmental
#' values from given observations at their year of sampling. Although the function
#' was developped mainly for CHELSA (https://chelsa-climate.org/), it can easily be
#' used for other type of data sources if following the same file format (i.e. single
#' raster files, informed environments and year in the URL or file names)
#'
#' @param env_layers Vector of file names or URL links of desired environmental layers. The
#' given strings must have the identity of the variable and the year informed.
#' @param env_id Vector. What are the unique names (as specified in 'env_layers') of the target
#' variables?
#' @param year_id Vector. What are the unique years (as specified in 'env_layers') we want the
#' environment from?
#' @param obs A data.frame with three columns (1-Longitude, 2-Latitude and 3-year of sampling).
#' The observations needs to be from the same CRS as the 'env_layers'.
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
#' # URL links
#' url.base = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/"
#' url.sprint1 = sprintf(paste0("tmax/CHELSAcruts_tmax_10_%d_V.1.0.tif"),1901:1911)
#' url.sprint2 = sprintf(paste0("tmin/CHELSAcruts_tmin_10_%d_V.1.0.tif"),1906:1921)
#' URL.links = paste0(url.base,c(url.sprint1,url.sprint2))
#'  
#' # Dummy observations
#' random.year = (1901:1921)[sample(1:21,nrow(xy.ppm),replace=TRUE)]
#' obs = data.frame(xy.ppm[,c("x","y")],year=random.year)
#' 
#' # Run function
#' sp.env = get_env(env_layers = URL.links,
#'                  env_id = c("tmax","tmin"),
#'                  year_id = 1901:1921,
#'                  obs = obs,
#'                  d.path = getwd())
#' 
#' @export
#' 
get_env = function(env_layers=NULL,
				   env_id=NULL,
				   year_id=NULL,
				   obs=data.frame(),
				   d.path=getwd())
{
	# URL or files
	if (grepl("https://",env_layers[1])) {

		# Main path to download and/or extract the rasters
		out.path = paste0(d.path,"/",gsub(".*/","",env_layers))

		# Loop
		lapply(1:length(env_layers), function(x){

			# Download
			ftp.file = try(download.file(env_layers[x],out.path[x],mode="wb"))
			if (is(ftp.file,"try-error"))
			{
				while(is(ftp.file,"try-error"))
				{
					cat("Reconnecting...","\n",sep="")
					ftp.file = try(download.file(env_layers[x],out.path,mode="wb"))
				}
			}
		})
	} else {
		out.path = paste0(d.path,"/",env_layers)
	}

	# Open the rasters
	our.ras = lapply(out.path, function(x) rast(x))

	# Prepare out data.frame
	out.df = as.data.frame(matrix(nrow=nrow(obs),ncol=ncol(obs)+length(env_id)))
	out.df[,1:3] = obs
	names(out.df) = c(names(obs),env_id)

	# Years that needs extraction
	ext.years = unique(out.df[,3])

	# Loop over environmental variable
	for (i in 1:length(env_id)) {

		# Get the right rasters
		list.p1 = gsub(".*/","",sapply(our.ras,sources))
		id.ras = our.ras[grepl(env_id[i],list.p1)]
		list.p2 = list.p1[grepl(env_id[i],list.p1)]

		# Loop over years
		for (j in 1:length(ext.years)) {

			# Observations and predictors from this year
			year.obs = out.df[out.df[,3]%in%ext.years[j],1:2]
			year.ras = id.ras[grepl(ext.years[j],list.p2)]
			if (length(year.ras)==0) {next} else {year.ras=year.ras[[1]]}

			# Extract infos
			env.info = extract(year.ras,year.obs)
			out.df[out.df[,3]%in%ext.years[j],3+i] = env.info[,!names(env.info)%in%"ID"]
		}
	}

	# Return
	return(out.df)
}

