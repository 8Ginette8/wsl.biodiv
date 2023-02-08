# wsl.biodiv

Developing a toolbox to efficiently handle SDMs in collaborative projects. Package aims at providing varied and flexible model fitting and tools for SDMs. Meta information should be stored (author, date, taxon,...) efficiently. Model evaluation and prediction should flexibly handle various fitted model objects. Evaluation metrics for presence-only and presence-absence models should be mutually implemented. Filtering and testing tools should help the user to find adequate predictors to apply meaningful models. Poisson Point Process models (PPPM) implementation should remain user-friendly. Regularization and variable selection should be implemented in this framework. Efficient methods to correct sampling bias in model fitting should be implements (pseudo-absences sempling strategies, environmental bias correction, bias covariate correction...)

# Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/wsl.biodiv")
library(wsl.biodiv)
```

# Example

## Ensemble model
### Data preparation

Load data:
``` r
# Take anguilla data set from dismo package
data("Anguilla_train")
vrs = c("SegSumT","USRainDays","USSlope")
env = Anguilla_train[,vrs]

# Alps data
data(AlpineConvention_lonlat)
data(exrst)
rst = rst[[1:6]]
data(xy_ppm)
mypoints = xy.ppm[,c("x","y")]
```

### Run ensemble

Custom parameters for the ensemble model:
``` r
# Formulas
form.glm = as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
form.gam = as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
form.gbm = as.formula(Presence ~ .)
feat = c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")

# All options
modinp = list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE,weight=TRUE),
   multi("gbm",list(formula=form.gbm,distribution = "bernoulli",interaction.depth = 1,shrinkage=.01,n.trees = 3500),"gbm-simple"), 
   multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE,weight=TRUE),
   multi("maxent",list(args=feat),"mxe-simple"),
   multi("randomForest",list(formula=form.gbm,ntree=500,maxnodes=NULL),"waud1"))
```

Calibrate ensemble model:
``` r
modi5 = wsl.flex(pa=Anguilla_train$Angaus,
               env_vars = env,
               taxon = "Angaus",
               replicatetype = "cv",
               reps = 3,
               strata = sample(1:3,nrow(env),replace=TRUE),
               project = "multitest",
               mod_args = modinp)
summary(modi5)
```

### Evaluate & predict

Evaluate and display results:
``` r
# Evaluate the model
eval5 = wsl.evaluate.pa(modi5,crit = "maxTSS")

# Get outputs or evaluation summary
eval5
summary(eval5)

# Get thresholds
thr.5 = get_thres(eval5, mean=FALSE)
```

Let's predict now (also works when 'predat' is a raster):
``` r
# Make some predictions (works also with Raster objects)
pred4 = wsl.predict.pa(modi5,predat=env)
pred5 = wsl.predict.pa(modi5,predat=env,thres=thr.5)
```

## Point process models (PPMs)
### Data preparation

Define a mask of your study area, to set a window and sample quadrature points:
``` r
# Define mask
maskR = mask(rst[[1]],shp.lonlat)

# Run 'wsl.ppm.window' function
wind = wsl.ppm.window(mask = maskR,
                      val = 1,
                      owin = TRUE)

# nDefine quadrature points for 'wsl.ppmGlasso'
quadG1 = wsl.quadrature(mask = maskR,
                        area.win = wind,
                        random = FALSE,
                        lasso = TRUE,
                        env_vars = rst)

# Define your environments
envG = raster::extract(rst,mypoints)
```

### Block cross-validation (BCV)

``` r
# Spatial block cross-validation
to_b_xy = rbind(mypoints,quadG1@coords)
toSamp = c(rep(1,nrow(mypoints)),rep(0,nrow(quadG1@coords)))
block_cv_xy = make_blocks(nstrat = 5, df = to_b_xy, nclusters = 10, pres = toSamp)
 
# Environmental block cross-validation
to_b_env = rbind(envG,quadG1@Qenv[,-1])
block_cv_env = make_blocks(nstrat = 5, df = to_b_env, nclusters = 10, pres = toSamp)
```

### Run PPMs

Fit a PPM with an Elastic Net regularization:
``` r
ppm.lasso = wsl.ppmGlasso(pres = mypoints,
                       quadPoints = quadG1,
                       asurface = raster::area(shp.lonlat)/1000,
                       env_vars = envG,
                       taxon = "species_eg1",
                       replicatetype = "cv",
                       reps = 5,
                       strata = NA,
                       save=FALSE,
                       project = "lasso_eg1",
                       path = NA,
                       poly = TRUE,
                       lasso = TRUE,
                       alpha = 0.5,             # 0.5 = Elastic Net here
                       type.measure = "mse",    # Other regularization parameters
                       standardize = TRUE,      # ....
                       nfolds = 5,              # ....
                       nlambda = 100)           # ....see cv.glmnet() for more details
summary(ppm.lasso)
```

Fit a simple PPM (without any regularization, nor polynomial terms) using environmental block cross-validation:
``` r
ppm.simple = wsl.ppmGlasso(pres = mypoints,
                       quadPoints = quadG1,
                       asurface = raster::area(shp.lonlat)/1000,
                       env_vars = envG,
                       taxon = "species_eg2",
                       replicatetype = "block-cv",
                       reps = 5,
                       strata = block_cv_env,
                       save = FALSE,
                       project = "lasso_eg2",
                       path = NA,
                       poly = FALSE,
                       lasso = FALSE)
summary(ppm.simple)
```

### Evaluate & predict

Evaluation example using Boyce index:
``` r
eval.lasso = wsl.evaluate.pres(x = ppm.lasso,
                               env_vars = rst,
                               speedup = TRUE)
summary(eval.lasso)
```

Evaluation example using other binary metrics:
``` r
eval.simple = wsl.evaluate.pa(x = ppm.simple,
                              crit = "maxTSS",
                              pres_only = TRUE)
summary(eval.simple)
```

Evaluation example by resetting a potential fitted bias covariate to 0 values (e.g. to remove sampling bias):
``` r
eval.bias = wsl.evaluate.pa(x = ppm.simple,
                            crit = "maxTSS",
                            pres_only = TRUE,
                            bias_cov = c(1,1,1,1,1,0))
summary(eval.bias)
```

Get calculated thresholds (mean may be chosen):
``` r
get_thres(eval.lasso, mean = FALSE)
get_thres(eval.simple, mean = TRUE)
get_thres(eval.bias, mean = FALSE)
```

Now we can predict:
``` r
# e.g. without using the thresholds --> species 'abundances'
pred.lasso = wsl.predict.pres(x = ppm.lasso,
                         predat = rst,
                         raster = TRUE)
                         
# e.g. using the thresholds --> species presences/absences
pred.simple = wsl.predict.pres(x = ppm.simple,
                         predat = rst,
                         thres = get_thres(eval.simple,mean=FALSE),
                         raster = TRUE)
                         
# e.g. or resetting a potential fitted bias covariate to 0 values:
pred.bias = wsl.predict.pres(x = ppm.simple,
                             predat = rst,
                             thres = get_thres(eval.bias,mean=FALSE),
                             raster = TRUE,
                             bias_cov = c(1,1,1,1,1,0))
```

Finally let's see what distributions we obtain across the European Alps by plotting:
``` r
par(mfrow=c(2,3))
sapply(1:5,function(x) plot(pred.lasso@predictions[[x]][[1]]))
```
![image](https://user-images.githubusercontent.com/43674773/217249539-61522c6b-3f8e-4779-9fa2-956549bb0397.png)

``` r
par(mfrow=c(2,3))
sapply(1:5,function(x) plot(pred.simple@predictions[[x]][[1]]))
```
![image](https://user-images.githubusercontent.com/43674773/217249813-94ba9338-6f74-495c-bc9e-b750e1f2399c.png)

``` r
par(mfrow=c(2,3))
sapply(1:5,function(x) plot(pred.bias@predictions[[x]][[1]]))
```
![image](https://user-images.githubusercontent.com/43674773/217249893-ee6c9c7c-9e26-4ffe-a1b7-f1debc69224b.png)

# Citations

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O., Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species distribution projections under climate change. Journal of Biogeography, 47(1), 130-142. doi: <a href="https://doi.org/10.1111/jbi.13734">10.1111/jbi.13734</a>

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>
