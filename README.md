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

Load package & data
``` r
# Package
library(wsl.biodiv)

# Take anguilla data set from dismo package
data("Anguilla_train")
vrs=c("SegSumT","USRainDays","USSlope")
env=Anguilla_train[,vrs]
```

Custom parameters for the ensemble model
``` r
# Formulas
form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
form.gbm=as.formula(Presence ~ .)
feat=c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")

# All options
modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE,weight=TRUE),
   multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE,weight=TRUE),
   multi("maxent",list(args=feat),"mxe-simple"),
   multi("randomForest",list(formula=form.gbm,ntree=500,maxnodes=NULL),"waud1"),
   multi("glm",list(formula=form.glm.2,family="binomial"),"glm-lin",step=TRUE,weight=TRUE))
```

Calibrate ensemble model
``` r
modi5=wsl.flex(pa=Anguilla_train$Angaus,
               env_vars = env,
               taxon="Angaus",
               replicatetype="block-cv",
               reps=3,
               strata=sample(1:3,nrow(env),replace=TRUE),
               project="multitest",
               mod_args=modinp)
```

Evaluate and display
``` r
# Evaluate the model
eval5<-wsl.evaluate.pa(modi5,crit="pp=op")

# Get outputs or evaluation summary
eval5
summary(eval5)
```

Let's predict now
``` r

```

## Point process model (PPM) lasso

``` r
...
```

# Citations

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O., Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species distribution projections under climate change. Journal of Biogeography, 47(1), 130-142. doi: <a href="https://doi.org/10.1111/jbi.13734">10.1111/jbi.13734</a>

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>
