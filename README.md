# wsl_biodiv

2018-05-08 V0.1

## Authors

Niklaus Zimmermann  
Rafael Wüest  
Philipp Brun (philipp.brun at wsl.ch)  

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

## Functions and classes

### Model fitting

#### Main functions
- wsl.glm
- wsl.gam
- wsl.maxent
- wsl.gbm

Flexibly fit various types of functions but let framework take care of resampling, meta-info storage, and file saving

- wsl.multi.fit

As above but a list with any models/parametrizations can be supplied for which predictions can be made.

- summary.wsl.fit

Print a summary of wsl.fit objects

#### Sub functions
- preps

Check input data, collect meta information, take care of data subsetting  

- multi

Generate a multi.input object for the wsl.multi.fit function (see below)

- hde

Avoid functions from printing unneccessary stuff

#### Classes

- wsl.fit

Generic class to store model output, testing data subsets, and meta information

- multi.input

Simple class to store model specifications for wsl.multi.fit function

### Model evaluation

### Model prediction
