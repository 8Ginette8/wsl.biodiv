### =========================================================================
### BOYCE INDEX
### =========================================================================
#' Boyce index model evaluation sub-function I
#'
#' Not to be called directly by the user
#'
#' @author 'ecospat.boyce' from the 'ecospat' package
#' @export
ecoBoyce=function (fit, obs, nclass = 0, window.w = "default", res = 100, 
          PEplot = TRUE) 
{
    if (class(fit) == "RasterLayer") {
        if (class(obs) == "data.frame" | class(obs) == "matrix") {
            obs <- extract(fit, obs)
        }
        fit <- getValues(fit)
        fit <- fit[!is.na(fit)]
    }
    if (window.w == "default") {
        window.w <- (max(fit) - min(fit))/10
    }
    interval <- c(min(fit), max(fit))
    mini <- interval[1]
    maxi <- interval[2]
    if (nclass == 0) {
        vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - 
                                                                    mini - window.w)/res)
        vec.mov[res + 1] <- vec.mov[res + 1] + 1
        interval <- cbind(vec.mov, vec.mov + window.w)
    }
    else if (length(nclass) > 1) {
        vec.mov <- c(mini, nclass)
        interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
    else if (nclass > 0 & length(nclass) < 2) {
        vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
    }
    f <- apply(interval, 1, boycei, obs, fit)
    to.keep <- which(f != "NaN")
    f <- f[to.keep]
    if (length(f) < 2) {
        b <- NA
    }
    else {
        r <- c(1:length(f))[f != c(f[-1], FALSE)]
        b <- cor(f[r], vec.mov[to.keep][r], method = "spearman")
    }
    HS <- apply(interval, 1, sum)/2
    HS[length(HS)] <- HS[length(HS)] - 1
    HS <- HS[to.keep]
    if (PEplot == TRUE) {
        plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", 
             col = "grey", cex = 0.75)
        points(HS[r], f[r], pch = 19, cex = 0.75)
    }
    results <- list(F.ratio = f, Spearman.cor = round(b, 3), 
                    HS = HS)
    return(results)
}
