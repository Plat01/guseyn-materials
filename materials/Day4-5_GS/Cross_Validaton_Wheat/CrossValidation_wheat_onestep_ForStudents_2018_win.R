###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                               Day 5 : Genomic Selection                                 ###
###                      -Case study #2: Prediction of yield in wheat: Cross-validation-    ###
###                                          ----                                           ###
###          R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, INP-ENSAT           ###
###                  amended from Dr Julio Isidro Sanchez & Dr Deniz Akdemir                ###
###                                                                                         ###
###############################################################################################

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
# use of constraint 'sum-to-zero' for ANOVA
options(contrasts=c('contr.sum','contr.poly'))

#############################################################################################################
## Loading of the R packages needed for the analysis.
library(rrBLUP)
library(lme4) 

###################
# Data import in R
###################

load("WheatData.yield.RData")

# Coding of the Amatfunc to compute the kinship matrix
Amatfunc<-function(MarkerMat){
  #MarkerMat is -1 0 1 coded
  pks<-colMeans(MarkerMat+1)/2
  W<-scale(MarkerMat, center=T, scale=F)
  c<-2*sum(pks*(1-pks))
  Amat<-tcrossprod(W)/c
  rownames(Amat)<-colnames(Amat)<-rownames(MarkerMat)
  return(Amat)
}
Amat.wheat<-Amatfunc(genodata)
Amat.wheat[1:5,1:5]
str(Amat.wheat)
dim(Amat.wheat) 


## preparing for cross-validation procedure
dim(genodata)
nlines <- nrow(genodata)

set.seed(1234)  ## initialize random number generator so that everyone in the audience will have the same 'random' numbers
nfolds <- 5 ## we decide to randomly split the dataset in 5 subsets


folds <- sample(1:nfolds, nlines, replace = TRUE)  ## assignation of each accession to 1 of the 5 subsets
#Binomial sampling: we take a line we assign it to a fold and we put it back. It can end up in another fold. We make a total of 337 samplings.
#On the other hand, some lines may not be sampled.
#Redo 50 times the cross validation.

folds

table(folds)  ## number of lines in each subset: 62+76+76+62+61=337

corrrBLUPvec<-c()  #prepare to put 5 correlations between obs. and predicted values into the vector

corGBLUPvec<-c()  #prepare to put 5 correlations into the vector

corBayesBvec<-c()

## -------------- the loop for cross-validation --------------

#Leave one out 5 times (5 folds = 5 subgroups)

for (rep in 1:1){  
  
  ## you may decide to run more than one cross-validation procedure

    ## in that case you need to sample accessions HERE and disable set.seed()

    for (i in 1:nfolds){

        print(paste("-------  this is  fold:", i, '  -------'))
        ## make sure that the order of factor is the same in the different datatsets (pheno and geno)
        phenodata.yield$Taxa <- factor(as.character(phenodata.yield$Taxa), levels = rownames(genodata))

        ## the lines that will belong to the test set in each loop
        samplelines<-rownames(genodata)[folds!=i]

        ## genotypic and phenotypic data of training set
        genodatatrain<-genodata[rownames(genodata) %in% samplelines,]
        phenodatatrain<-phenodata.yield[phenodata.yield$Taxa%in%samplelines,]

        
        ## genotypic and phenotypic (not usual.. because the aim is to predict them) data of test set
        genodatatest<-genodata[!(rownames(genodata)%in%samplelines),]
        phenodatatest<-phenodata.yield[!(phenodata.yield$Taxa%in%samplelines),]

        
        length(levels(phenodatatrain$Taxa))
        
        Xtrain<-model.matrix(~-1+Env,data=phenodatatrain)
        Xtest<-model.matrix(~-1+Env,data=phenodatatest)
        
        Ztrain<-model.matrix(~-1+Taxa,data=phenodatatrain)
        Ztest<-model.matrix(~-1+Taxa,data=phenodatatest)
        
        dim(Ztrain)
        dim(Ztest)

####rrBLUP
##        library(rrBLUP)
        print(paste("computing RR-BLUP model for fold", i))
        
        model.onestep.rrBLUP <- mixed.solve( y = phenodatatrain$yield,
                                            X = Xtrain,
                                            Z = Ztrain %*% genodata,
                                            K = NULL)   ## based on additive values of markers : RR-BLUP
        
        markereffects.onestep.rrBLUP <- matrix(model.onestep.rrBLUP$u, ncol=1)

        
########After you get marker effects we can get GEBVs
        GEBVstrain.RRBLUP.onestep <- genodatatrain %*% markereffects.onestep.rrBLUP  ## use mol. markers as random effects

        GEBVstest.RRBLUP.onestep <- genodatatest %*% markereffects.onestep.rrBLUP

######

        GEBVs.RRBLUP.onestep <- as.matrix(genodata) %*% markereffects.onestep.rrBLUP

########
        corrrBLUPvec <- c(corrrBLUPvec,cor(phenodatatest$yield, Ztest %*% GEBVs.RRBLUP.onestep))
        print(paste("Corr. RR-BLUP, fold", i, ":", cor(phenodatatest$yield, Ztest %*% GEBVs.RRBLUP.onestep)))

        
#GBLUP avec K=Amat.wheat
        print(paste("computing G-BLUP model for fold", i))
        
        model.GBLUP.onestep <- mixed.solve(y = phenodatatrain$yield,
                                           X = Xtrain,
                                           Z = Ztrain,
                                           K = Amat.wheat)  ## using realized kinship : G-BLUP

        
        GEBVstrain.gblup.onestep <- Ztrain %*% model.GBLUP.onestep$u  ## use Ztrain (accessions)  as random effects

        GEBVstest.gblup.onestep <- Ztest %*% model.GBLUP.onestep$u

######
        
        GEBVs.gblup.onestep <- model.GBLUP.onestep$u

########
        corGBLUPvec <- c(corGBLUPvec, cor(phenodatatest$yield, Ztest %*% GEBVs.gblup.onestep))
        print(paste("Corr. G-BLUP, fold", i, ":", cor(phenodatatest$yield, Ztest %*% GEBVs.gblup.onestep)))
        
        ### BAYESIAN LASSO ONE STEP
        library(BGLR)
        ?BGLR
        model.BayesB.onestep<-BGLR(phenodatatrain$yield,
                                   ETA=list(list(X=model.matrix(~-1+Env,data=phenodatatrain), model="FIXED"), 
                                            list(X=Ztrain%*%genodata,model="BayesB")),nIter=5000, burnIn=2000
                                   ,verbose=FALSE)
        summary(model.BayesB.onestep)
        
        marker.effects.BayesB.onestep<- model.BayesB.onestep$ETA[[2]]$b
        ########After you get marker effects we can get GEBVs
        GEBVstrain.BayesB.onestep<-genodatatrain%*%marker.effects.BayesB.onestep
        
        
        GEBVstest.BayesB.onestep<-genodatatest%*%marker.effects.BayesB.onestep
        
        ######
        
        GEBVs.BayesB.onestep<-as.matrix(genodata)%*%marker.effects.BayesB.onestep
        
        ########
        corBayesBvec<-c(corBayesBvec,cor(phenodatatest$yield,Ztest%*%GEBVs.BayesB.onestep))
        
        }

}


x11()
boxplot(cbind(corrrBLUPvec, corGBLUPvec, corBayesBvec))#To pick the best (i.e. more accurate) model
