###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                                 Day 4: Genomic Selection                                ###
###                                                                                         ###
###  Case study #1: Prediction of yield depending on Environments(BLUE)and genotypes(BLUP)- ###
###  Dataset created from experimental design described in Roumet et al., 2012, New Phytol  ###
###                                          ----                                           ###
###        R Script edited by Pr Laurent Gentzbittel & C�cile Ben, INP-ENSAT                ###
###############################################################################################

######### CASE STUDY PRESENTATION ##########################################################
###Prediction of genetic value of sunflower cultivars for yield based on the phenotypic results of a Multi-Environment Trial (MET)
# 6 sunflower genotypes -"NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"- 
# in 5 environments near Toulouse (France) -"Odars", "Verfeil", "Villefranche", "Fourquevaux"- 
# evaluated in a VERY UNBALANCED design.
# during 3 years -2003, 2005, 2007 -


######### PREPARATION OF THE WORKING INTERFACE IN R ##########################################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space:clean workspace, clean memory 
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# garbage collector. To use to free memory
gc() 
#############################################################################################################

## Loading of the R packages needed for the analysis.
library(rrBLUP)
library(lme4) 

###################
# Data import in R
###################
#########
## EX. 1: Computation of BLUEs & BLUPs to predict Yield using environnement and Genotype WITHOUT MARKERS
#########

sunflower <- read.table("GS1_SunflowerWeed.csv", sep=";", header=TRUE, dec=".")
sunflower #Dataset 'inspired' from Roumet et al., New Phytol, 2012
str(sunflower)
#6 genotypes of cultivated sunflowers, 5 environments, 3 years, 
#6 genotypes = markers in real GS datasets

sunflower$YEAR<-as.factor(sunflower$YEAR)
str(sunflower)

##Vector of observations------------------------------------------------------------------------
y <-sunflower$YLD
y #Dim (12,1)

## Design matrix restricted to environnements------------------------------------------------------------------------

#Let's create X matrix by hand.
# we use the constraint 'set-first-to-zero'

Xmanual <- matrix(c(1,1,0,0,0,   #µ + Effect of experimental field (i.e Odars) for Measure1 
              1,1,0,0,0,	 #µ + Effect of experimental field (i.e Odars) for Measure2 
              1,0,1,0,0,	 #µ + Effect of experimental field (i.e. Verfeil) for Measure3 
              1,0,1,0,0,   #µ + Effect of experimental field (i.e. Verfeil) for Measure4 
              1,0,0,1,0,   #µ + Effect of experimental field (i.e. Villefranche)for Measure5 
              1,0,0,1,0,
              1,0,0,1,0,	
              1,0,0,0,1,
              1,0,0,0,1,
              1,0,0,0,1,
              1,0,0,0,0,
              1,0,0,0,0),
            byrow=TRUE,nrow=12,ncol=5)
dimnames(Xmanual)<-list(c("Measure1","Measure2", "Measure3", "Measure4", "Measure5", "Measure6", "Measure7", "Measure8", "Measure9", "Measure10", "Measure11", "Measure12"),
                    c("Intercept", "Odars", "Verfeil", "Villefranche", "Fourquevaux"))
Xmanual
#another way to create X is 
X <- model.matrix(~sunflower$ENV)
X_env_year <- model.matrix(~sunflower$ENV+sunflower$YEAR)
#Rq: In this dataset, no blocks, no replicates by enviro and year so no possibility to take
# into account the Site x Year interaction.
#If it was the case: X_envXyear <- model.matrix(~sunflower$ENV*sunflower$YEAR)


#2 solutions to solve linear syst
#Sol. 1: ANOVA Constraint : sum to 0: Sum of alpha effects (alpha Env1 + alpha Env2 + alpha Env3=0) but difficult to program
#solution 2: ANOVA Constraint : set to 0: effect of the 1st modalité by alphabetical order (here, Balma) set to 0.
# Effect of this enviro = intercept.
#The X matrix created with the function contains  (Env-1) columns: Balma column is not present.
# Balma is the first alphabetical Env.
X
#Identical to the X matrix we designed manually except that by default the function model.matrix order
# Factor level by alphabetical order and remove the 1st level

EnvBetas=solve((t(X)%*%X))%*%t(X)%*%y
EnvBetas
#Effect of Balma = 0 due to constraint "set-to-zero"
#Estimated mean = 1479.50000
Env_yearBetas=solve((t(X_env_year)%*%X_env_year))%*%t(X_env_year)%*%y
Env_yearBetas
#Beta_Balma equal 0 and Beta_2003 equal 0
#Balma Estimated mean = intercept and 2003 Estimated mean = intercept.
#If Site X Year interaction in the model then intercept = Estimated mean for Balma_2003.

#Decomposition of the solving of this linear system
#solve(t(X)%*%X)#Calcule l'inverse de (t(X)%*%X)
#solve(t(X)%*%X)%*%(t(X)%*%X)#  Identity Matrix identité "à l'erreur près" proving what is computed by (t(X)%*%X)
#t(X)%*%y

## Mixed Linear Model with 'Genotypes' as random effects:
#Creation of Z matrix of random effects and of µ vector of BLUP
sunflower
Zmanual<- matrix(c(1,0,0,0,0,0,   #Ind.1=µ1
                   0,1,0,0,0,0,  	#Ind.2=µ2
                   1,0,0,0,0,0,		#Ind.3=µ3
                   0,1,0,0,0,0,
                   0,1,0,0,0,0,
                   0,0,1,0,0,0, 
                   0,0,0,1,0,0,
                   1,0,0,0,0,0,	
                   0,0,0,0,1,0,
                   0,0,0,0,0,1,
                   0,1,0,0,0,0,
                   0,1,0,0,0,0
),byrow=TRUE,nrow=12,ncol=6)
dimnames(Zmanual)<-list(c("Measure1","Measure2", "Measure3", "Measure4", "Measure5", "Measure6", "Measure7", "Measure8", "Measure9", "Measure10", "Measure11", "Measure12"),
                  c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"))

Zmanual

## ------------------------------------------------------------------------
Z <- model.matrix(~sunflower$Geno-1)
#By default, model.matrix add an intercept. However, for random effect, no dependencies so no intercept.  
#By default, Genotypes in alphabetical order in Z matrix 
Z

## Calculation of BLUP from centered values of phenotypes (y -mean(y))
# without taking into account environment effects
# here, we want to fit the following model : Prédiction model for BLUP: y2 = mean(y) + Zu + E 
#' Solving the system for random effects needs to take into account
#' the VarCov matrix of the random effects (i.e. that will be later the kinship matrix). 
#' In this first analysis, we consider that there is no kinship between genotypes and we will use identity matrix. I = diag(6)

y2 <- (sunflower$YLD) - mean(sunflower$YLD) #Obs. Mean = 1502.333
uhatSimple <- solve(t(Z) %*% Z + 1*diag(6)) %*% t(Z) %*% y2
#1*diag(6)=lambdaG with sigmaE/sigmaU arbitrarily set to 1

uhatSimple


#' # Fitting a M.E.T. with a mixed model.

#' MultiEnvironment Trials: Estimation of genetic values but NOT of the breeding values (additive)
# ------------------------------------------------------------------------
# Mixed model: Fixed effect for Environnements & random effects for genotypes

library(lme4)
modellmer<-lmer(sunflower$YLD ~ sunflower$ENV + (1|sunflower$Geno))
#By default, lmer uses G = SigmaU*I & R = SigmaE*I 
modellmer 
#Compare the BLUE for environnements from this model and the ones computed with fixed effects only:
EnvBetas

ranef(modellmer) 
#Compare the BLUP for genotypes from this model and the ones computed with random effects only:
uhatSimple

#To visualize the X design matrix and Z random effect matrix
getME(modellmer,'X')
t(as.matrix(getME(modellmer,'Zt')))#unexpected complex syntax...Sorry...

#In "modellmer", Please see:

#Groups     Name        Std.Dev.
# sunflower$Geno (Intercept) 54.36   
# Residual               17.96  
#VarG=54.36^2 VarE=17.96^2 h²= Broad Sense heritability
# except if genotypes are pure lines and design is balanced (narrow-sense heritability).
## Caution: If unbalanced design as here, the environnement effects and genotype effects are not
## independent, so this does not provide an accurate Broad sense heritability.


##estimate of broad-sense heritability = square of accuracy
as.data.frame(VarCorr(modellmer))#vcov for heritability calculation
# (i.e. VarG & Var E) are in column 4
(varianceestimates<-as.data.frame(VarCorr(modellmer))[,4])
h2hat<-varianceestimates[1]/sum(varianceestimates)
h2hat
varianceestimates[1]/varianceestimates[2] #SigmaE/SigmaU 9.160897


#' # Estimating genetic effects, taking into account relationships between varieties
#' 


###
#Reminder: In the previous example, we do not know anything about the relationships between the varieties
#So we used the identity matrix I in the model.

#' ## Fitting a random model using pedigrees.
## Model P-BLUP (P for pedigree) based on a (fake) pedigree matrix

Amat<- matrix(c(1,0.5,0.25,0,0.25,0,
                0.5,1,0,0.25,0.5,0.25,
                0.25,0,1,0.5,0.5,0,
                0,0.25,0.5,1,0,0.25,
                0.25,0.5,0.5,0,1,0.5,
                0,0.25,0,0.25,0.5,1
               ),byrow=TRUE,nrow=6,ncol=6)
dimnames(Amat)<-list(c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"),
                        c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"))
Amat

#Clone: 1; Full sibs: 0.5; Half sibs: 0.25; 1st cousins: 1/8; Parent-offsprings: 0.5; GrandParents/Grandchild: 0.25; Aunt(Uncle)-Nieces(nephew): 0.25
breedingvalue.withA <- solve(t(Z) %*% Z + 1*solve(Amat)) %*% t(Z) %*% y2
breedingvalue.withA
#Compare the BLUP for genotypes from this model and the ones computed with random effects only.
uhatSimple
# Compare the BLUP for genotypes from this model and the ones computed with a Mixed model: Fixed effect for Environnements & random effects for genotypes
ranef(modellmer) 

#' ## Fitting a random  model using a kinship matrix.
#'

#' ### Marker genotypes and imputation of missing genotypic data.
#' 
#########
## Taking into account of a Kinship matrix
#FAKE DATA: 6 genotypes and 10 markers - Imputation of genotypes and calculation of the Kinship matrix
#########

MarkersMatrix = read.table("MToImp.csv", sep=";", header=TRUE) #Reads your file
dim(MarkersMatrix)  #Tell you the dimensions of object MarkersMatrix in rows and columns
MarkersMatrix
MForImp= as.matrix(MarkersMatrix) #Creates matrix from imported data
is.matrix(MForImp)
dimnames(MForImp)<-list(c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"), 
                        c("M1","M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
                        )
dim(MForImp)
MForImp #Coding for alleles depends on their frequency :
# minor hmz allele : +1, htz: 0, major hmz allele: -1

library(rrBLUP)
Gimputed <- A.mat(MForImp, max.missing = NULL,
                 impute.method="mean", return.imputed = T)
#For each marker,imputation of the missing data by the mean of th observed genotypes. 
Gimputed$imputed #To print the matrix of imputed genotypes.
# Ex: marker 2: 3*1/4=0.75; marker3 (2*-1+1)4=-0.25

#Another possible method to impute the missing genotype: EM algorithm
Gimputed.EM <- A.mat(MForImp, max.missing = NULL,
                    impute.method="EM", return.imputed = T)
Gimputed.EM$imputed

# rrBLUP::Amat also computes a kinship matrix.
# Caution : This kinship matrix in $A can present some negative values.
# It corresponds to 'first G' estimator of VanRaden (2008)
# negative genomic relationships : This	implies that individuals i and j are
# more distinct than an average pair of individuals in the data
# Fixing negative estimates of relationships to	0 is **wrong** praxis.


#' ### Computing the kinship matrix ``by hand''.
#'

#Canonical calculation of the Kinship matrix
p <- colSums((Gimputed$imputed + 1)/ 12)  # 6 diploid accessions thus a maximum of 12 alleles
p2 <- 2 * sum(p*(1-p))
Gmatrix <- crossprod(t(Gimputed$imputed))/p2
# this matrix is usally different form that computed with rrBLUP::Amat


#' ### Prediction of breeding values.
#'
#'

# this is a random-only model (w/o fixed environmental effects)
breedingvalue.withG <- solve(t(Z) %*% Z + 1*solve(Gmatrix)) %*% t(Z) %*% y2  # by hand
breedingvalue.withG
# compare with uhatSimple (same model BUT identity matrix as G) or to
# breedingvalue.withA (same model BUT (fake) pedigree matrix as G)


# this is the mixed model (fixed environments effects and random genetic values)
#library(rrBLUP)
GBLUPmodel <- mixed.solve(y = y,X = X,Z = Z, K = Gmatrix, method = "REML")  # using rrBLUP::mixed.solve
#The function mixed.solve also estimate sigmaU and sigmaE
#G=K*SigmaU and R=I*SigmaE
GBLUPmodel
GBLUPmodel$Vu / (GBLUPmodel$Vu + GBLUPmodel$Ve)

# Restricted maximum likelihood (REML) because the G matrix is unknown but estimated.
# For manual solving of the Henderson's equation, we use least squares method
# which only apply if we exactly know the G matrix (assumption emitted in this case)
# with rrBLUP impossible to code for least squares method
# LL : Lod likelihood value : when chosing between 2 models, take the one with the highest LL.


#' # Predicting phenotypes, only knowing genetic relationships.

#'
#' Prediction of the phenotype of a new line (which was not phenotyped)
#' on the basis of its genotypic kinship with varieties that have been evaluated in the field
## ------------------------------------------------------------------------
###

set.seed(1)
# generating genotypes for three new lines
newline1 <- sample(c(-1,0,1),10,replace=T)
newline2 <- sample(c(-1,0,1),10,replace=T)
newline3 <- sample(c(-1,0,1),10,replace=T)
Mnew <- rbind(Gimputed$imputed,newline1,newline2,newline3)
rownames(Mnew) <- c("Atomic",
                    "Santiago",
                    "NK Countri",
                    "Boogy",
                    "Aurasol",
                    "PR64H41",
                    "newline1",
                    "newline2",
                    "newline3")
Mnew

## Calculation of the kinship matrix: relationships of the new lines with the first set of accessions
# ------------------------------------------------------------------------
Mnew+1
p <-colSums(Mnew+1)/ 18  # now 9 diploid genotypes, so 18 alleles
sum(2*p*(1-p))
Gnew <- Mnew %*% t(Mnew)/ (sum(2*p*(1-p)))   ## Van Raden's 'first G'
Gnew
dim(Gnew) #  pay attention to the dimensions.

# ------------------------------------------------------------------------
(sunflower2 <- sunflower)

sunflower2$Geno <- factor(as.character(sunflower$Geno),
                          levels=rownames(Gnew))
Znew <- model.matrix(~ sunflower2$Geno-1)
dim(Znew)
## THIS dimension IS VERY important : 12 realized observations and 9 genetic effects to estimate
# 6 genotypes with phenotypes and 3 genotypes W/O phenotypes.


## Prediction using a random-only model: Caution ! we do not take into account the effect of the environnements (no X matrix)
# ------------------------------------------------------------------------
breedingvalue.new <- solve(t(Znew) %*% Znew + 1 * solve(Gnew)) %*% t(Znew) %*% y2
# here : SigmaE/SigmaU = 1 by hypothesis
breedingvalue.new 


## prediction using the full mixed model.
# ------------------------------------------------------------------------
mixed.solve(y = y, X = X, Z = Znew, K = Gnew, method = "REML")  # mixed.solve estimates SigmaU and SigmaE
