# wsl_biodiv

2018-05-17 V0.1

## Authors

Philipp Brun (philipp.brun at wsl.ch)  
Rafael Wüest  
Niklaus Zimmermann  

Dynamic Macroecology Group  
Swiss Federal Research Institute WSL  

## Aim

Developing a toolbox to efficiently handle SDMs in collaborative projects.

- Model fitting should be kept as flexible as possible
- Meta information should be stored (author, date, taxon,...)
- Files should be saved in a structured way
- Relationships between saved objects should be preserved
- Model evaluation should flexibly handle various fitted model objects
- Model prediction should flexibly handle various fitted model objects

## Main functions

### Fit models

- __wsl.glm__
- __wsl.gam__
- __wsl.maxent__
- __wsl.gbm__
- __wsl.flex__

Flexibly fit various types of functions but let framework take care of resampling, meta-info storage, and file saving. wsl.flex basically allows supplying any possible model, however, there may be problems with prediction/evaluation for exotic functions.

#### Arguments

- pa: vector with presence/absence values
- env_vars: data.frame with environmental predictors
- taxon: name of the taxon for which models are fitted
- replicatetype: (how) should replicates be generated? may be 'none', 'splitsample', 'cv' 'block-cv'
- reps: number of replicates
- strata: a numeric vector of the same length as observations with integers separating cross validation replicates (used when replicatetype='block-cv')
- save: should the model be saved in a structured way? (not implemented yet)
- project: character indicating the name of the project within which the models are run (later used to define saving directories)
- path: where to save? (not implemented yet)
- step: (for glms and gams only) should the models be updated with the step function?
- mod_tag: (not in wsl.flex) label for the current model
- mod_args: list with elements of class 'multi.input' which specify models to be fitted in wsl.flex

#### Output
Object of class wsl.fit including slots for meta info, testing data for evaluation, and model objects

### Summary of model fits

- __summary.wsl.fit__

Print summary of wsl.fit objects

### Evaluate models

- __wsl.evaluate__

Assess several model skill metrics for all models in a wsl.fit object. Currently AUC, RMSE, TSS, PPV, Accuracy, and Cohen's Kappa are evaluated. Furthermore, the threshold applied is returned.

#### Arguments

- x: a wsl.fit object
- tester: data.frame with testing data (only mandatory if replicatetype='none' was chosen when models were fitted)
- threshold: vector of the same length as number of models chosen with custom thresholds for model evaluation. for wsl.flex outputs the thresholds have to be labelled with the same names provided to models
- crit: which threshold criterion should be considered? Currently 'pp=op' (predicted prevalence = observed prevalence), 'maxTSS' (threshold yielding maximum TSS), and 'external' (thresholds manually supplied) are possible 

#### Output
Object of class wsl.evaluation with slots for meta info, and model performance estimates 

### Summary of model evaluation

- __summary.wsl.evaluation__

Print summary of wsl.evaluation objects

### Obtain thresholds from evaluation object

- __get_thres__

Extracts thresholds from wsl.evaluation objects and names them so they can be fed to the wsl.predict function. At the moment only averages over replicates can be obtained.

#### Arguments

- x: a wsl.evaluation object

#### Output

A named vector with one value per model.

### Make predictions

- __wsl.predict__

Make predictions with all models from a wsl.fit object. If thresholds are supplied binary predictions are made, otherwise continuous predictions are returned.

#### Arguments

- x: a wsl.fit object
- predat: data.frame with points for which predictions should be made
- thres:  vector of the same length as number of models chosen with custom thresholds for model evaluation. for wsl.flex outputs the thresholds have to be labelled with the same names provided to models

#### Output
Object of class wsl.prediction with slots for meta info, and model predictions 

### Define model settings for wsl.flex function

- __multi__

Create a multi.input object that efficiently stores model specifications

#### Arguments

- mod: a character with the name of the function to be called. E.g. "gam"
- args: a list with arguments to be passed to the function specified in mod
- tag: character with name for model set-up
- step: should step function be applied to update model

#### Output
Object of class multi

### Block-wise split data into training and testing

- __make_blocks__

Creates a stratum vector based on a data.frame with n columns. If the data.frame has one column strata are created based on clusters separated by quantiles. If the data.frame has two or more columns, strata ere created based on k-medoid clusters (function 'pam' from package cluster). Instead of a data.frame also the argument 'npoints' can be provided, then groups are created by random sampling. An opitimization algorithm (function 'gridSearch' from package NMOF) optimizes for equal stratum sizes. 

#### Arguments

- nstrata: number of approximately equal-sized classes to separate groups in block-cross validation
- df: data.frame with n columns containing critera for cluster building. Not necessary if argument npoints is supplied
- nclusters: number of clusters based on which strata should be built. Minimum the same number as starta, maxuimum nrow(df)/10
- npoints: optional argument if 'df' is not supplied. For how many points should random sampling be made?

#### Output
Vector of length nrow(df) or npoints, with integers representing different strata

### Sample pseudo-absences proportional to presence point distribution

- __prop.sampling__

Uses the 'density' function from the package spatstat to create a density surface of the supplied point pattern and samples pseudo-absences from this density-distribution.

#### Arguments

- points: matrix or data.frame with column names 'x' and 'y' assumed to be on the same scale and metric distances (m, km,...)
- nsamples: number of pseudoabsences to be generated
- res: resolution of the denstiy grid from which pseudo-absences are drawn. default (1) corresponds to 1000 cells on the x-axis
- ...: arguments passed on to the 'density' function from package spatstat. Note in particluar the argument 'adjust' which controls the kernel size in the density interpolation.

#### Output
nsamples x 2 matrix with drawn psuedo-absences

## Helper functions (not to be called directly by the user)

- __preps__

Check input data, collect meta information, take care of data subsetting. Called by model fitting functions.

- __hde__

Avoid functions from printing unneccessary stuff. Called by model fitting functions.

- __ceval__

Do the actual model evaluations

- __preva.meta__

Generate meta information for prediction and evaluation

- __prd__

Correctly feed the predict functions depending on model type (glm, gbm, maxent...)

- __optme__

Optimization function to create equal-sized strata in the 'make_blocks' function
