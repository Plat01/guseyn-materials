###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                                 Day 4: Genomic Selection                                ###
###                                                                                         ###
###                      -Case study #2: Prediction of yield in wheat                       ###
###                 depending on Environments(BLUE)and SNP effects(BLUP)-                   ###
###                                          ----                                           ###
###          R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, INP-ENSAT           ###
###                  amended from Dr Julio Isidro Sanchez & Dr Deniz Akdemir                ###
###                                                                                         ###
###############################################################################################

# Preparing working space and data.

## libraries.
######### PREPARATION OF THE WORKING INTERFACE IN R ##########################################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed,
# no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space:clean workspace, clean memory 
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# garbage collector. To use to free memory
gc() 
# use of constraint 'sum-to-zero' for ANOVA
options(contrasts=c('contr.sum','contr.poly'))

#############################################################################################################

## Loading of the R packages needed for the analysis.
library(rrBLUP)
library(lme4)

## Datasets : genotypes for training set and test set, and phenotypes for training set.

###################
# Data import in R
###################

genodata <- as.matrix(read.table("./genodata.csv", sep = ";", header = TRUE))

#The ``genodata`` dataframe contains the genotypic data (3735 SNPs) for ALL the 337 lines: training set (250 lines) and test set (87 lines).
# line's names are indicated as rownames.
dim(genodata)  # lines as rows, genotypes as columns
## We will use the lines names on this file as the reference
## to keep the different matrix consistent.

######
###Phenotypic and genotypic data for training population
######

phenodatatrain <- read.table("./phenodatatrain.csv", sep=";", header=T)
dim(phenodatatrain)
phenodatatrain$Taxa <- factor(as.character(phenodatatrain$Taxa),
                              levels = rownames(genodata))  ## <-- key syntax (see sunflower example)
# Making sure that the levels of pheno.train are the same than the original geno.train
# you must be sure that Z matrix row and col are in the same order row.names and col.names
# of the kinship matrix 

str(phenodatatrain)          ## 337 levels for Taxa (training + test sets)
unique(phenodatatrain$Taxa)  ## but only 250 lines in the training dataset
phenodatatrain
# Yield measurements for a total population of 250 lines in 4 environments
# (combinations of 2012_Davis/2013_Imperial X Dry/Irr)
# Some lines are tested in all environments (e.g. H0800310),
# some only in two of them (e.g. SELKIRK in 2012_Davis X Dry/Irr)
#IN TOTAL, 694 measures

genodatatrain <- read.table("genodatatrain.csv", sep=";", header=T)
dim(genodatatrain)  # lines as rows, genotypes as columns
#str(genodatatrain) #Training population contain 250 individuals that have been genotyped for 3735 markers
# subset of genodata for the training set.



######
###Genotypic data for test population
######

#Test population contain 87 individuals that have been genotyped for 3735 markers
genodatatest <- read.table("./genodatatest.csv", sep=";", header=T)
str(genodatatest)

dim(genodatatrain)  # lines in rows, SNPs in columns
dim(genodatatest)   # lines in rows, SNPs in columns
#Train+Test=337

# some data organisation that we will need latter
GenodataTest <- genodatatest
GenodataTest$Taxa <- rownames(genodatatest)
GenodataTest$Taxa <- factor(as.character(GenodataTest$Taxa), levels = rownames(genodata)) 


# Building the model matrices.

#We choose to work with matrices without the constant term ($\mu$) so that there is no intercept.

######
###X and Z matrices for training population
######

# For fixed effect. We remove the intercept.
Xtrain <- model.matrix(~ -1 + Env, data = phenodatatrain)
dim(Xtrain)  ## is it expected ?
head(Xtrain)
head(phenodatatrain)


Ztrain <- model.matrix(~ -1 + Taxa, data = phenodatatrain)
str(Ztrain)
dim(Ztrain) #Measure in rows and Taxa in columns
#337 columns (although only  250 lines in the training set) *BUT* the columns of Ztrain corresponding to test set are NULL. 
#As an example, lines ``UC1710`` and``MTHW1069`` are in the test set (see ``rownames(genodatatest)``).
sum(Ztrain[, "TaxaUC1710"]) ## 0 as expected
sum(Ztrain[, "TaxaMTHW1069"])   ## 0 as expected

# head(Ztrain)  # if you need to have a look.

######
### Prediction with a G-BLUP model and a kinship matrix
######

#Use markers to model genetic relatedness among lines then use relatedness to make predictions.
## Computation of the kinship matrix

# Coding of the Amatfunc to compute the kinship matrix
# it corresponds to VanRaden's 2008
Amatfunc<-function(MarkerMat){
  #MarkerMat is -1 0 1 coded
  pks <- colMeans(MarkerMat+1)/2
  W <- scale(MarkerMat, center=T, scale=F)
  c <- 2*sum(pks*(1-pks))
  Amat <- tcrossprod(W)/c
  rownames(Amat) <- colnames(Amat) <- rownames(MarkerMat)
  return(Amat)
}

#we use ALL genotypic data to compute the kinship matrix because we compute the marker-based genetic relationships between all the lines in the training set and all the lines in the test set. This is the basis of the VarCov of (random) genetic effects : $G = \sigma^2_u K$

Amat.wheat <- Amatfunc(genodata)
Amat.wheat[1:5,1:5]
str(Amat.wheat)
dim(Amat.wheat)  ## check if these are the appropriate dimensions


## Computation of GEBVs using the G-BLUP.

#This model ***is*** **the prediction model**.
# --------------------------------------------------------------------------------
# Ztrain has 337 columns, with the columns for test set set to zero. 
# the training set and the test set are 'connected' by the K matrix.
# THIS ALLOWS FOR THE BLUP OF the TEST SET TO BE COMPUTED
model.GBLUP.onestep <- mixed.solve(y = phenodatatrain$yield,
                                   X = Xtrain,
                                   Z = Ztrain,
                                   K = Amat.wheat)
#
# --------------------------------------------------------------------------------
model.GBLUP.onestep

# Notice that the u vector of random effects has the right dimension.
GEBVs.gblup.onestep <- model.GBLUP.onestep$u
dim(GEBVs.gblup.onestep)  ## is it expected ?


# BLUP for a couple of lines from the training set :
model.GBLUP.onestep$u["GLADIUS"]
model.GBLUP.onestep$u["CARBERRY"]

# BLUP for a couple of lines from the test set :
model.GBLUP.onestep$u["UC1710"]
model.GBLUP.onestep$u["MTHW1069"]


# If needed, Prediction with a G-BLUP model including an identity matrix:
# K set to NULL means 'K = I'
model.GBLUP.onestep.Kid <- mixed.solve(y = phenodatatrain$yield,
                                       X = Xtrain,
                                       Z = Ztrain) ## <-- K is not specified
model.GBLUP.onestep.Kid


# Getting the Estimation of Genomic Estimated Breeding Values (GEBVs) for lines within the training population
# this is simply the matrix product. Because Ztrain has 1/0 for training set and 0 elsewhere for test set
# it gives the BLUP of the training set
GEBVstrain.gblup.onestep <- Ztrain %*% model.GBLUP.onestep$u
dim(GEBVstrain.gblup.onestep)
GEBVstrain.gblup.onestep
dim(Ztrain)


#Getting the Prediction of Genomic Estimated Breeding Values (GEBVs) for lines within the test population
Ztest <- model.matrix(~ -1 + Taxa, data = GenodataTest)
dim(Ztest)
## In total 337 columns, but columns corresponding to the training set are NULL :
sum(Ztest[, "TaxaGLADIUS"])
sum(Ztest[, "TaxaCARBERRY"])


# this is simply the matrix product. Because Ztest has 1/0 for test set and 0 elsewhere for training set
GEBVstest.gblup.onestep <- Ztest %*% model.GBLUP.onestep$u
dim(GEBVstest.gblup.onestep)
GEBVstest.gblup.onestep




## Calculation of the accuracy of the prediction model. We will compare EBV's with observed yields.

#In fact, we had phenotypic data (... :-) )

phenodatatest <- read.table("./phenodatatest.csv", sep=";", header=T)
dim(phenodatatest)
phenodatatest$Taxa <- factor(as.character(phenodatatest$Taxa),
                             levels = rownames(genodata))  # <-- key syntax 
#Making sure that the levels of pheno.train are the same than the original geno.
str(phenodatatest)
unique(phenodatatest$Taxa)  ## 87 lines in the data set but 337 lines at all (training + test sets)


Ztest2 <- model.matrix(~ -1 + Taxa, data = phenodatatest)
dim(Ztest2)


## correlations btween observed and predicted yield:
cor(phenodatatest$yield, Ztest2 %*% GEBVs.gblup.onestep)  ## No GxE in the model
model.GBLUP.onestep$Vu
model.GBLUP.onestep$Ve
# but we have to account for heritability of the trait.
her.est <- (model.GBLUP.onestep$Vu/(model.GBLUP.onestep$Vu+model.GBLUP.onestep$Ve))
(accuracy_GBLUP_onestep<-cor(phenodatatest$yield, Ztest2 %*% GEBVs.gblup.onestep)/ sqrt(her.est))  #  <-- divide by heritability sqrt

######
### Genomic Selection with the RR-BLUP method.
######
#Use markers as predictors of genetic values then use estimated markers effects to make predictions.
#A model of genomic selection using Ridge-Regression uses the genotypes at markers as random effects -- and not the **u** for each genotype --. 
#the assumption is that there is no VarCov between the markers. K = I.

model.onestep.rrBLUP <- mixed.solve(y = phenodatatrain$yield,
                                    X = Xtrain,
                                    Z = Ztrain %*% genodata,  ## <-- genotypes
                                    K = NULL)
#Caution: K = NULL corresponds to K = I
 
# the dimension of Z matrix is much higher: there are 3735 random effects to be computed (instead of 337)
dim(Ztrain %*% genodata)
# please note that the genotypes of the test set are not used, because
# the columns of Ztrain corresponding to the test set are NULL (see above)

str(model.onestep.rrBLUP)

# effects of each of the 3735 SNPs
# need to load the position of the markers onto the genome
mrkMap <- read.table("./mrkMap.csv", sep = ";")

x11()
markereffects.onestep.rrBLUP<-matrix(model.onestep.rrBLUP$u, ncol=1)
plot(abs(markereffects.onestep.rrBLUP), ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=mrkMap$chrom,main='Marker Effects')



########After you get marker effects we can get GEBVs
GEBVstrain.RRBLUP.onestep <- as.matrix(genodatatrain) %*% markereffects.onestep.rrBLUP
# this is simply the sum of the effects of each allele at each of the 3735 SNPs for the training set
str(GEBVstrain.RRBLUP.onestep)
GEBVstrain.RRBLUP.onestep


GEBVstest.RRBLUP.onestep <- as.matrix(genodatatest) %*% markereffects.onestep.rrBLUP
# this is simply the sum of the effects of each allele at each of the 3735SNPs for the test set
str(GEBVstest.RRBLUP.onestep)
GEBVstest.RRBLUP.onestep
######

GEBVs.RRBLUP.onestep <- as.matrix(genodata) %*% markereffects.onestep.rrBLUP
# all GEBVs
str(GEBVs.RRBLUP.onestep)
GEBVs.RRBLUP.onestep
########


cor(phenodatatest$yield, Ztest2 %*% GEBVs.RRBLUP.onestep)

model.onestep.rrBLUP$Vu
model.onestep.rrBLUP$Ve

(her.est <- (model.onestep.rrBLUP$Vu/(model.onestep.rrBLUP$Vu+model.onestep.rrBLUP$Ve/3735)))
(accuracy_RRBLUP_onestep<-cor(phenodatatest$yield, Ztest2%*%GEBVs.RRBLUP.onestep)/sqrt(her.est))


####CORRELATION MATRIX

cor.mat<-cor(data.frame(rrBLUP=GEBVs.RRBLUP.onestep,
                        GBLUP=GEBVs.gblup.onestep))
cor.mat  ## conclusions ?

# all GEBVs, expected to be gaussian random effects
GEBVs.RRBLUP.onestep[order(GEBVs.RRBLUP.onestep)]
x11()
plot(GEBVs.RRBLUP.onestep[order(GEBVs.RRBLUP.onestep)])


# all training and test GEBVs, expected to be gaussian random effects
GEBVstest.RRBLUP.onestep[order(GEBVstest.RRBLUP.onestep)]

x11()
plot(GEBVstest.RRBLUP.onestep[order(GEBVstest.RRBLUP.onestep)])

x11()
plot(GEBVstrain.RRBLUP.onestep[order(GEBVstrain.RRBLUP.onestep)])

######
### BAYESIAN LASSO ONE STEP
######
##Use markers as predictors of genetic values then use estimated markers effects to make predictions.

library(BGLR)
?BGLR
model.BayesB.onestep<-BGLR(phenodatatrain$yield,
                           ETA=list(list(X=model.matrix(~-1+Env,data=phenodatatrain), model="FIXED"), 
                                    list(X=Ztrain%*%genodata,model="BayesB")),nIter=5000, burnIn=2000
                           ,verbose=FALSE)
#ETA: This is a two-level list used to specify the regression function (or linear predictor). 

summary(model.BayesB.onestep)

marker.effects.BayesB.onestep<- model.BayesB.onestep$ETA[[2]]$b

plot(abs(marker.effects.BayesB.onestep), ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=mrkMap$chrom,main='Marker Effects')

########After you get marker effects we can get GEBVs
GEBVstrain.BayesB.onestep<-as.matrix(genodatatrain)%*%marker.effects.BayesB.onestep


GEBVstest.BayesB.onestep<-as.matrix(genodatatest)%*%marker.effects.BayesB.onestep

######

GEBVs.BayesB.onestep<-as.matrix(genodata)%*%marker.effects.BayesB.onestep

########
cor(phenodatatest$yield,Ztest2%*%GEBVs.BayesB.onestep)

#PB!!!! No possible to compute accuracy!!! because model.BayesB.onestep$Vu does not exist!
(her.est <- (model.BayesB.onestep$Vu/(model.BayesB.onestep$Vu+model.BayesB.onestep$Ve/3735)))
(accuracy_BayesB_onestep<-cor(phenodatatest$yield, Ztest2%*%GEBVs.BayesB.onestep)/sqrt(her.est))

######
###RKHS
######

#Use markers to model genetic relatedness among lines then use relatedness to make predictions.

Eucl.distance<-as.matrix(dist(genodata))
h<-1/ncol(genodata)
Amat.RKHS<-exp(-h*(Eucl.distance)^2/(2))
Amat.RKHS[1:5,1:5]
Amat.wheat[1:5,1:5]

model.RKHS<-mixed.solve(y= phenodatatrain$yield, X=Xtrain, Z=Ztrain, K=Amat.RKHS)

GEBVstrain.RKHS<-Ztrain%*%model.RKHS$u
dim(GEBVstrain.RKHS)
GEBVstrain.RKHS

GEBVstest.RKHS<-Ztest%*%model.RKHS$u
dim(GEBVstest.RKHS)
GEBVstest.RKHS

GEBVs.RKHS<-model.RKHS$u

#Accuracy
(her.est <- (model.RKHS$Vu/(model.RKHS$Vu+model.RKHS$Ve/3735)))
(accuracy_RKHS_onestep<-cor(phenodatatest$yield, Ztest2%*%GEBVs.RKHS)/sqrt(her.est))

### Picking h based on likelihood statistics: Pb Visiblement pas de diff. avec le modèle précédent!
nmarkers<-ncol(genodata)
hvec<-rep(1/nmarkers,5)*c(1/3,1/2,1,2,3)
hvec

RHKS.hfun<-function(h){
  Amat.RKHS<-exp(-h*(Eucl.distance)^2/(2))
  model.RKHS<-mixed.solve(y= phenodatatrain$yield, X=Xtrain, Z=Ztrain, K=Amat.RKHS)
  return(model.RKHS$LL)
}

RKHS.models<-sapply(hvec,RHKS.hfun,simplify=TRUE)
RKHS.models
h.best<-hvec[which.max(RKHS.models)]
plot(RKHS.models)

#Use this h now to see if the correlations are better after we selected.

Amat.RKHS.best.h<-exp(-h.best*(Eucl.distance)^2/(2))
Amat.RKHS.best.h[1:5,1:5]


model.RKHS.best.h.onestep<-mixed.solve(y= phenodatatrain$yield, X=Xtrain, Z=Ztrain, 
                                       K=Amat.RKHS.best.h)

GEBVstrain.RKHS.best.h.onestep<-Ztrain%*%model.RKHS.best.h.onestep$u


GEBVstest.RKHS.best.h.onestep<-Ztest%*%model.RKHS.best.h.onestep$u

######

GEBVs.RKHS.best.h.onestep<-model.RKHS.best.h.onestep$u

########
#cor(phenodatatest$yield, Ztest2%*%GEBVs.RKHS.best.h.onestep)
(her.est <- (model.RKHS.best.h.onestep$Vu/(model.RKHS.best.h.onestep$Vu+model.RKHS.best.h.onestep$Ve/3735)))
(accuracy_RKHS.best.h_onestep<-cor(phenodatatest$yield, Ztest2%*%GEBVs.RKHS.best.h.onestep)/sqrt(her.est))

####CORRELATION MATRIX

cor.mat<-cor(data.frame(GBLUP=GEBVstest.gblup.onestep,
                        rrBLUP=GEBVstest.RRBLUP.onestep,
                        BayesB=GEBVstest.BayesB.onestep,
                        RKHS=GEBVstest.RKHS.best.h.onestep))
cor.mat
