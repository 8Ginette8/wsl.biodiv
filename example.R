# ###########################################################################
# Get an overview over the data and figure out which traits make most sense
# to be used
#
# $Date: 2017-11-29
#
# Authors: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 
# ###########################################################################

### =========================================================================
### Initialise system
### =========================================================================
# Start with house cleaning
rm(list = ls()); graphics.off()

### =========================================================================
### Get ID
### =========================================================================

setwd("R:/brunp/cluster_jobs/SDM_correlation/Src/Tools/wsl_biodiv/")

### =========================================================================
### Initialise system
### =========================================================================

# Load necessary Packages
library(mgcv)
library(dismo)
library(MASS)
library(gbm)
library(randomForest)
library(ROCR)
library(dismo)
library(cluster)
library(NMOF)
library(spatstat)

# source functions
scr=list.files("functions/",full.names=TRUE)

invisible(lapply(scr, source))

### =========================================================================
### Prepare data
### =========================================================================

# Take anguilla data set from dismo package
data("Anguilla_train")

# Run wsl glm
vrs=c("SegSumT","USRainDays","USSlope")

env=Anguilla_train[,vrs]

### =========================================================================
### Check out wsl.glm
### =========================================================================

form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))

# Try out wsl.glm funcion
modi1=wsl.glm(pa=Anguilla_train$Angaus,
              env_vars = env,
              taxon="Angaus",
              replicatetype="cv",
              reps=5,
              project="prototest",
              mod_tag="test-glm",
              formula=form.glm,
              family="binomial",
              step=TRUE)

# Try out custom summary function
summary(modi1)

# Access glm object of first replicate
summary(modi1@fits$replicate_01$`test-glm`)

# Evaluate the model
eval1=wsl.evaluate(modi1)

# Get evaluation summary
summary(eval1)

### =========================================================================
### Check out wsl.gam
### =========================================================================

form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))

# Try out wsl.glm funcion
modi2=wsl.gam(pa=Anguilla_train$Angaus,
             env_vars = env,
             taxon="Angaus",
             replicatetype="splitsample",
             reps=3,
             project="prototest",
             mod_tag="test-gam",
             formula=form.gam,
             family="binomial",
             step=FALSE)

# Try out custom summary function
summary(modi2)

# Access glm object of first replicate
summary(modi2@fits$replicate_01$`test-gam`)

# Evaluate the model
eval2=wsl.evaluate(modi2,crit="max")

# Get evaluation summary
summary(eval2)

### =========================================================================
### Check out wsl.gbm
### =========================================================================

form.gbm=as.formula(Presence ~ .)

# Try out wsl.glm funcion
modi3=wsl.gbm(pa=Anguilla_train$Angaus,
              env_vars = env,
              taxon="Angaus",
              replicatetype="none",
              reps=1,
              project="prototest",
              mod_tag="test-brt",
              formula= form.gbm,
              distribution = "bernoulli",
              interaction.depth = 1,
              shrinkage=.01,
              n.trees = 3500)

# Try out custom summary function
summary(modi3)

# Access glm object of first replicate
summary(modi3@fits$replicate_01$`test-brt`)

# Prepare external testing data
tste=data.frame(Presence=Anguilla_train$Angaus,env)

# Evaluate the model
eval3=wsl.evaluate(modi3,crit="max",tester=tste)

# Get evaluation summary
summary(eval3)

### =========================================================================
### Check out wsl.maxent
### =========================================================================

feat=c("linear=true","quadratic=true","hinge=true","product=true","threshold=false")

# Try out wsl.glm funcion
modi4=wsl.maxent(pa=Anguilla_train$Angaus,
                 env_vars = env,
                 taxon="Angaus",
                 replicatetype="block-cv",
                 reps=3,
                 strata=sample(1:3,nrow(env),replace=TRUE),
                 project="prototest",
                 mod_tag="test-mxe",
                 args=feat)

# Try out custom summary function
summary(modi4)

# Access glm object of first replicate
summary(modi4@fits$replicate_01$`test-mxe`)

# Define external threshold
thmxe=c(`test-mxe`=0.5)

# Evaluate the model
eval4=wsl.evaluate(modi4,crit="external",thres=thmxe)

# Get evaluation summary
summary(eval4)

### =========================================================================
### Check out wsl.multi.fit
### =========================================================================

form.glm.2=as.formula(paste("Presence~",paste(vrs,collapse="+")))

modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE),
            multi("gbm",list(formula=form.gbm,
                             distribution = "bernoulli",
                             interaction.depth = 1,
                             shrinkage=.01,
                             n.trees = 3500),"gbm-simple"),
            multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE),
            multi("maxent",list(args=feat),"mxe-simple"),
            multi("randomForest",list(formula=form.gbm,ntree=500),"waud"),
            multi("glm",list(formula=form.glm.2,family="binomial"),"glm-lin",step=TRUE))

# Try out wsl.glm funcion
modi5=wsl.multi.fit(pa=Anguilla_train$Angaus,
                    env_vars = env,
                    taxon="Angaus",
                    replicatetype="block-cv",
                    reps=3,
                    strata=sample(1:3,nrow(env),replace=TRUE),
                    project="multitest",
                    mod_args=modinp)

# Try out custom summary function
summary(modi5)

# Access glm object of first replicate
summary(modi5@fits$replicate_01$`glm-simple`)

# Evaluate the model
eval5<-wsl.evaluate(modi5,crit="max")

# Get evaluation summary
summary(eval5)

### =========================================================================
### Make some predictions
### =========================================================================

# Make prediction
pred4=wsl.predict(modi4,predat=env,thres=c(0.5))

pred5=wsl.predict(modi5,predat=env)

### =========================================================================
### Test out block generation function
### =========================================================================

# No layers supplied
strt.1=make_blocks(npoints=1000)
table(strt.1)

# Stratified by 1d layer a
strt.2=make_blocks(df=env[,1,drop=F],nstrat=5,nclusters=5)
table(strt.2)

# Stratified by 1d layer b
strt.3=make_blocks(df=env[,1,drop=F],nstrat=5,nclusters=15)
table(strt.3)

# Stratified by 2d layer a
strt.4=make_blocks(df=env[,c(1,3)],nstrat=3,nclusters=3)
table(strt.4)

# Stratified by 2d layer b
strt.5=make_blocks(df=env[,c(1,3)],nstrat=5,nclusters=50)
table(strt.5)

# Stratified by 3d layer
strt.6=make_blocks(df=env[,1:3],nstrat=5,nclusters=50)
table(strt.6)

par(mfrow=c(3,2))
plot(env[,c(1,3)],col=strt.1)
plot(env[,c(1,3)],col=strt.2)
plot(env[,c(1,3)],col=strt.3)
plot(env[,c(1,3)],col=strt.4)
plot(env[,c(1,3)],col=strt.5)
plot(env[,c(1,3)],col=strt.6)

### =========================================================================
### Test proportional sampling
### =========================================================================

load("mypts.RData")

# Samples with decreasing density kernel length
pseu.1=prop.sampling(points=my.pts,nsamples=10000,adj=1,res=1)
pseu.2=prop.sampling(points=my.pts,nsamples=10000,adj=.5,res=1)
pseu.3=prop.sampling(points=my.pts,nsamples=10000,adj=.2,res=1)
pseu.4=prop.sampling(points=my.pts,nsamples=10000,adj=.1,res=1)

par(mfrow=c(2,2))

plot(pseu.1,pch=16,col="red",cex=.5)
points(my.pts)
plot(pseu.2,pch=16,col="red",cex=.5)
points(my.pts)
plot(pseu.3,pch=16,col="red",cex=.5)
points(my.pts)
plot(pseu.4,pch=16,col="red",cex=.5)
points(my.pts)

