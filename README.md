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

- wsl.glm
- wsl.gam
- wsl.maxent
- wsl.gbm
- wsl.multi.fit

Flexibly fit various types of functions but let framework take care of resampling, meta-info storage, and file saving. wsl.multi.fit basically allows supplying any possible model, however, there may be problems with prediction/evaluation for exotic functions.

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
- mod_tag: (not in wsl.multi.fit) label for the current model
- mod_args: list with elements of class 'multi.input' which specify models to be fitted in wsl.multi.fit

#### Output
Object of class wsl.fit including slots for meta info, testing data for evaluation, and model objects

### Summary of model fits

- summary.wsl.fit

Print summary of wsl.fit objects

### Evaluate models

- wsl.evaluate

Assess several model skill metrics for all models in a wsl.fit object. Currently AUC, RMSE, TSS, PPV, Accuracy, and Cohen's Kappa are evaluated.

#### Arguments

- x: a wsl.fit object
- tester: data.frame with testing data (only mandatory if replicatetype='none' was chosen when models were fitted)
- threshold: vector of the same length as number of models chosen with custom thresholds for model evaluation. for wsl.multi.fit outputs the thresholds have to be labelled with the same names provided to models
- crit: which threshold criterion should be considered? Currently 'pp=op' (predicted prevalence = observed prevalence), 'max' maximum value of among all cutoffs, and 'external' (thresholds manually supplied) are possible 

#### Output
Object of class wsl.evaluation with slots for meta info, and model performance estimates 

### Summary of model evaluation

- summary.wsl.evaluation

Print summary of wsl.evaluation objects

### Make predictions

- wsl.predict

Make predictions with all models from a wsl.fit object. If thresholds are supplied binary predictions are made, otherwise continuous predictions are returned.

#### Arguments

- x: a wsl.fit object
- predat: data.frame with points for which predictions should be made
- thres:  vector of the same length as number of models chosen with custom thresholds for model evaluation. for wsl.multi.fit outputs the thresholds have to be labelled with the same names provided to models

#### Output
Object of class wsl.prediction with slots for meta info, and model predictions 

### Define model settings for wsl.multi.fit function

- multi
Create a multi.input object that efficiently stores model specifications

#### Arguments

- mod: a character with the name of the function to be called. E.g. "gam"
- args: a list with arguments to be passed to the function specified in mod
- tag: character with name for model set-up
- step: should step function be applied to update model

#### Output
Object of class multi

#### Helper functions (not to be called directly by the user)
- preps
Check input data, collect meta information, take care of data subsetting. Called by model fitting functions.

- hde
Avoid functions from printing unneccessary stuff. Called by model fitting functions.

- ceval
Do the actual model evaluations

- preva.meta
Generate meta information for prediction and evaluation

- prd
Correctly feed the predict functions depending on model type (glm, gbm, maxent...)
