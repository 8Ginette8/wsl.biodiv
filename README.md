# wsl_biodiv

Developing a toolbox to efficiently handle SDMs in collaborative projects. Package aims at providing varied and flexible model fitting and tools for SDMs. Meta information should be stored (author, date, taxon,...) efficiently. Model evaluation and prediction should flexibly handle various fitted model objects. Evaluation metrics for presence-only and presence-absence models should be mutually implemented. Filtering and testing tools should help the user to find adequate predictors to apply meaningful models. Poisson Point Process models (PPPM) implementation should remain user-friendly. Regularization and variable selection should be implemented in this framework. Efficient methods to correct sampling bias in model fitting should be implements (pseudo-absences sempling strategies, environmental bias correction, bias covariate correction...)

# Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/wsl.biodiv")
library(wsl.biodiv)
```

# Example

## Ensemble model

Load package & data:
``` r
# Package
library(wsl.biodiv)

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

Custom parameters for the ensemble model:
``` r
# Formulas
form.glm = as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
form.gam = as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
form.gbm = as.formula(Presence ~ .)
feat = c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")

# All options
modinp = list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE,weight=TRUE),
   multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE,weight=TRUE),
   multi("maxent",list(args=feat),"mxe-simple"),
   multi("randomForest",list(formula=form.gbm,ntree=500,maxnodes=NULL),"waud1"))
```

Calibrate ensemble model:
``` r
modi5 = wsl.flex(pa=Anguilla_train$Angaus,
               env_vars = env,
               taxon="Angaus",
               replicatetype="block-cv",
               reps=3,
               strata=sample(1:3,nrow(env),replace=TRUE),
               project="multitest",
               mod_args=modinp)
summary(modi5)
```

Evaluate and display results:
``` r
# Evaluate the model
eval5 = wsl.evaluate.pa(modi5,crit="maxTSS")

# Get outputs or evaluation summary
eval5
summary(eval5)
```

Let's predict now:
``` r
# Make some predictions (works also with Raster objects)
pred4=wsl.predict.pa(modi4,predat=env)
pred5=wsl.predict.pa(modi5,predat=env,thres=thr.5)
```

## Point process models (PPM)
### Data preparation

Define a mask of your study area, to set a window and sample quadrature points:
``` r
# Define mask
maskR = mask(rst[[1]],shp.lonlat)

# Run 'wsl.ppm.window' function
wind = wsl.ppm.window(mask = maskR,
                      val = 1,
                      owin = TRUE)

# Define random  quadrature points for 'wsl.ppmGlasso'
quadG1 = wsl.quadrature(mask = maskR,
                        area.win = wind,
                        random = FALSE,
                        lasso = TRUE,
                        env_vars = rst)
```

### Block cross-validation (BCV)

Spatial BCV:
``` r
to_b_xy = rbind(mypoints,quadG1@coords)
toSamp = c(rep(1,nrow(mypoints)),rep(0,nrow(quadG1@coords)))
block_cv_xy = make_blocks(nstrat = 5, df = to_b_xy, nclusters = 10, pres = toSamp)
```

Environmental BCV:
``` r
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

Fit a simple PPM (without any regularization) using block cross-validation:
``` r
ppm.simple = wsl.ppmGlasso(pres = mypoints,
                       quadPoints = quadG1,
                       asurface = raster::area(shp.lonlat)/1000,
                       env_vars = envG,
                       taxon = "species_eg2",
                       replicatetype = "cv",
                       reps = 5,
                       strata = NA,
                       save = FALSE,
                       project = "lasso_eg2",
                       path = NA,
                       poly = FALSE,
                       lasso = FALSE)
summary(ppm.simple)
```

# Citations

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O., Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species distribution projections under climate change. Journal of Biogeography, 47(1), 130-142. doi: <a href="https://doi.org/10.1111/jbi.13734">10.1111/jbi.13734</a>

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>
