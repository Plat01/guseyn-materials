###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                                Day 1: Augmented Block Design                            ###
###                                                                                         ###
###           Case Study #4: Evaluation trials for new varieties registration               ### 
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, INP-ENSAT          ###
###############################################################################################

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'sum-to-zero' for ANOVAs
options(contrasts=c('contr.sum','contr.poly'))
#contr.sum: sum of effects of factors = 0

#############################################################################################################

## Classically, analysis is done with a fixed model
## Preliminary analysis with MLM is proposed in a second part


## Loading of the R packages needed for the analysis.

##library(car)    # Levene & de Brown test
library(agricolae)   # Newman-Keuls & de Tukey test
library(ggplot2)
library(tidyverse)

library(lmerTest)
library(lsmeans)   # For mean estimation



###################
# Data import in R
###################

Produc <- read.table("04_AugmentedBlockDesign_example1.csv",sep=";",dec=".",header=TRUE)
str(Produc)


### Experimental design?
with(Produc,
     table(Accession,Block))



######################
###################### Analysis with fixed model, all Accessions
######################

## type I SS

Ana1a <- aov(weight ~ Block + Accession , data = Produc)  # Accessions 'corrected' from Block
summary(Ana1a)


Ana1b <- aov(weight ~ Accession + Block , data = Produc)  # Blocks 'corrected' from Accessions  
summary(Ana1b)
## The two ANOVA tables are different because the design is unbalanced :
## some genotypes are absent in some Blocks 



###########
###################### Analysis on checks only

## fixed model, Accessions 'corrected' from Block
AnaChecks <- aov(weight ~ Block + Accession , data = subset(Produc,subset=Nature == "check"))
summary(AnaChecks)


###########
###################### Values of Accessions 'corrected from Block effects' using checks information over the design.

########### SS type I method

## Effects of Blocks, to correct the values of strains
model.tables(AnaChecks)
model.tables(AnaChecks)$tables$Block

### calculation of the data 'corrected from Block effects'
Produc$EffetBlock <-sapply(Produc$Block,switch,
                          B1=model.tables(AnaChecks)$tables$Block[1],
                          B2=model.tables(AnaChecks)$tables$Block[2],
                          B3=model.tables(AnaChecks)$tables$Block[3]
                          )

Produc$weightCorr <- Produc$weight - Produc$EffetBlock


## To check :
## ATTENTION : it works ONLY for checks, because the experimental design is balanced only for checks


Balanced1 <- aov(weight ~ Block + Accession, data = subset(Produc,subset=Nature == "check"))
summary(Balanced1)
Balanced2 <- aov(weight ~ Accession + Block, data = subset(Produc,subset=Nature == "check"))
summary(Balanced2)

Balanced3 <- aov(weightCorr ~ Block + Accession, data = subset(Produc,subset=Nature == "check"))
summary(Balanced3)
Balanced4 <- aov(weightCorr ~ Accession + Block, data = subset(Produc,subset=Nature == "check"))
summary(Balanced4)

# note Block SS in Balanced1/2 vs Balanced 3/4
# note Accession SS in Balanced1/2 vs Balanced 3/4
# Is it expected ?


###########  maeans adjusted for Block effect, using lsmeans :
lsmeans(Ana1a, ~ Accession)


## compare weightCorr et lsmeans values
## to understand WHAT are lsmeans :
aggregate(Produc$weightCorr ~ Produc$Accession, FUN='mean')




###########
###################### Analysis of the variability between checks, between checks and new Accessions, between new Accessions
## Technical (and laborious) part but important to get it.

# contrast:  checks vs new varieties
cm.checks <- c(-2,-2,-2,-2,1,1,1,1,1,1,1,1)
# contrasts:  among checks
cm.C1vsC2 <- c(1,-1,0,0,0,0,0,0,0,0,0,0)
cm.C1C2vsC3 <- c(1,1,-2,0,0,0,0,0,0,0,0,0)
cm.C1C2C3vsC4 <- c(1,1,1,-3,0,0,0,0,0,0,0,0)
# contrasts:  among new varieties
cm.V1vsV2 <- c(0,0,0,0,1,-1,0,0,0,0,0,0)
cm.V1V2vsV3 <- c(0,0,0,0,1,1,-2,0,0,0,0,0)
cm.V1V2V3vsV4 <- c(0,0,0,0,1,1,1,-3,0,0,0,0)
cm.V1V2V3V4vsV5 <- c(0,0,0,0,1,1,1,1,-4,0,0,0)
cm.V1V2V3V4V5vsV6 <- c(0,0,0,0,1,1,1,1,1,-5,0,0)
cm.V1V2V3V4V5V6vsV7 <- c(0,0,0,0,1,1,1,1,1,1,-6,0)
cm.V1V2V3V4V5V6V7vsV8 <- c(0,0,0,0,1,1,1,1,1,1,1,-7)

# Create a matrix that contains the comparisons : 
mat <- cbind(cm.checks,
             cm.C1vsC2,cm.C1C2vsC3,cm.C1C2C3vsC4,
             cm.V1vsV2,cm.V1V2vsV3,cm.V1V2V3vsV4,cm.V1V2V3V4vsV5,cm.V1V2V3V4V5vsV6,cm.V1V2V3V4V5V6vsV7,cm.V1V2V3V4V5V6V7vsV8
             )

# Apply that matrix to the data ("sum-to-zero")
contrasts(Produc$Accession) <- mat
contrasts(Produc$Accession)


# Compute ANOVA with new contrasts and show SSQ decompositions

## treatment (Accessions) adjusted:
AnaFixe1 <- aov( weight ~ Block + Accession, data = Produc)
## Block Adjusted
AnaFixe2 <- aov( weight ~ Accession + Block, data = Produc)

summary(AnaFixe1, split=list( Accession=list( "among checks" = c(2,3,4),
                                          "among varieties" = c(5,6,7,8,9,10,11),
                                          "checks vs varieties "= 1))
        )


summary(AnaFixe2, split=list( Accession=list( "among checks"=c(2,3,4),
                                          "among varieties"=c(5,6,7,8,9,10,11),
                                          "checks vs varieties "=1))
        )


# The summary command includes the option: split.
# The split option provides a list of factors where contrasts are stored.
# Within each factor, we also provide a list which includes the names of the contrasts (i.e. each column) stored in the contrasts matrix columns (i.e. mat).



###########
################### Conclusion : Which are the Accessions different from checks ?

c <- 4
r <- 3
ddl <- (r-1)*(c-1)
MSE <- summary(Ana1a)[[1]]$'Mean Sq'[3]  # 3rd term of the column MS from Ana1a summary

# SE of the difference between 2 adj means of selections in different Blocks (cf slides)
Sv <- sqrt(2*(c+1)*MSE/c)

#vSE of the difference between adjusted selection mean and check. This us the useful quantity to know !
Svc <- sqrt(((r+1)*(c+1)*MSE)/(r*c)) 

# We want to know the Accessions which are better than checks
# t threshold value : 
qt(0.95,ddl)

LSI <- qt(0.95,ddl) * Svc
LSI



###################
######################################  Analysis with agricolae
###################

# To be secure, we specify again constraints on ANOVA effects
options(contrasts=c("contr.sum", "contr.poly"))
## 
contrasts(Produc$Accession) <- "contr.sum"

with(Produc,
by(Accession, Block, as.character)
)


modelDAU1 <- DAU.test(Produc$Block, Produc$Accession, Produc$weight, method = "lsd", console = TRUE)
## The table of "ANOVA, Treatment Adjusted" is in Ana1a ( Block + Accessions)
## The table of "ANOVA, Block Adjusted" is in Ana1b (Accessions + Block)

### Please see the contrasts decomposition done in the models AnaFixe1 & AnaFixe2 to understand how agricolae works.

## Small numerical diff on Svc between agricolae et and by-hand computation

options(digits = 2)
modelDAU1$means
### Please see lsmeansand corrected means to xompare with values from agricolae.


## Detail : If we want 2-by-2 pairs comparison and not by group .
modelDAU2 <- DAU.test(Produc$Block, Produc$Accession, Produc$weight, method="lsd", group = FALSE, console = TRUE)
head(modelDAU2$comparison,8)

## In conclusion : agricolae uses a fixed model for ANOVA.
## We can do better than agricolae.






###########
############################################  Analysis in Mixed Linear Model (MLM), all Accessions
###########
# TO be secure, we specify again the constraints of ANOVA

contrasts(Produc$Accession) <- "contr.sum"

########### in this example, the number of Blocks is too small to
########### accuratly estimate their random effect. BUT we will pursue with this example
########### to compare with fixed models.


########### Case 1 : Blocks are random effects, genotypes are fixed effects

Ana2a <- lme4::lmer(weight ~ Accession + (1 | Block), data = Produc)
summary(Ana2a)
anova(Ana2a)     ## note that ddl of fixed effect are NOT specified (for some theoratical reasons)

Ana2b <- lmerTest::lmer(weight ~ Accession + (1 | Block), data = Produc)
summary(Ana2b)
anova(Ana2b)    ## note that ddl of fixed effect are different from those of the fixed model 
                ## see ana1a ou ana1b (fixed model)

## SCE is 'Block ajusted' *BUT* df are different so the power is higher because critical F are lower
qf(c(0.95, 0.99), 11, 6)  ## for fixed model
qf(c(0.95, 0.99), 11, 8)  ## for mixed model


## Analysis on checks only, in MLM
AnaChecks2 <- lmerTest::lmer(weight ~ Accession + (1 | Block), data = subset(Produc,subset = Nature == "check"))
summary(AnaChecks2)


##
# Xb 
Fix.BLUEs <- getME(AnaChecks2,'X') %*% fixef(AnaChecks2)
Fix.BLUEs
# Zu
Ran.BLUPs <- t(as.matrix(getME(AnaChecks2,'Zt'))) %*% unlist(ranef(AnaChecks2))  
Ran.BLUPs
# Xb + Zu
Predict <- Fix.BLUEs + Ran.BLUPs
Predict

## compare with result of predict()function of  model, rounding number to 6digits
all(round((Predict), 6) == round(predict(AnaChecks2),6))  ## comparison with estimates given by the mixed model. Fully logical !
all(round((Predict), 6) == round(predict(AnaChecks), 6))  ## comparison with estimates given by the fixed model. weird , no ?!?


lsmeans(Ana2a, ~ Accession)
fixef(Ana2a) + fixef(Ana2a)[1]
## for the fun : it lacks a value. How to compute it ?


AdjustMoys <- lsmeans::lsmeans(Ana2a ,  ## lmerTest is masking standard lsmeans. Unwanted !
                     pairwise ~ Accession, 
                     adjust = "tukey")

cld(AdjustMoys,
    alpha = 0.05 ,
    Letters = letters ,
    adjust = "tukey"
    )




########### Case 2 : Blocks are random effects, genotypes are random effects
## both random effects are uncorrelated


Ana3 <- lmerTest::lmer(weight ~ 1 + (1 | Block) + (1|Accession), data = Produc)
summary(Ana3)

anova(Ana3)    ## quite stupid : there is only the grand mean as fixed effect


## numerical values of the random effects : Accession and Block
ranef(Ana3)

## BLUP of the Accessions : We got them !!!
ranef(Ana3)$Accession + fixef(Ana3)




## prediction interval (not confidence interval !) of the BLUPs :

library(merTools)  ## warning need an object from lme4::lmer. Interesting library to explore

Ana3b <- lme4::lmer(weight ~ 1 + (1 | Block) + (1|Accession), data = Produc)
Randoms <- REsim(Ana3b, n.sims = 500)
## It provides the name of the grouping factor, the level of the factor we are obtaining an estimate for, the term in the model, and the mean, median, and standard deviation of the simulated values. We can use this to generate a caterpillar plot 

x11()
plotREsim(Randoms,
          labs = TRUE)
## the values that have a confidence interval that does not overlap zero are highlighted in black.
## You can modify the width of the interval by using the level parameter to plotREsim making wider or narrower confidence intervals based on your needs. 




## theory of the MLM :
## to print the different matrices of the model
extract.lmerModDesign(Ana3)#Please, execute first the function "extract.lmerModDesign" below.



##############################################################################################
##############################################################################################
##
## functions
##
##############################################################################################


##
##############################################################################################
### extract.lmerModDesign function should be defined in RLRsim, but did not work (20171104)

#' @importFrom stats model.response
#' @importFrom lme4 getME VarCorr
#' @rawNamespace 
#' if(getRversion() >= "3.3.0") {
#'   importFrom("stats", sigma)
#' } else {
#'   importFrom("lme4", sigma)
#' }
extract.lmerModDesign <- function(m) {
  X <- getME(m,"X")
  Z <- as.matrix(getME(m,"Z"))
  v <- VarCorr(m)
  resvar <- sigma(m)^2
  Sigma.l <- lapply(v,function(x) x/resvar) #Cov(b)/ Var(Error)
  k <- getME(m,"n_rtrms") #how many grouping factors
  q <- lapply(Sigma.l,NROW) #how many variance components in each grouping factor
  ## OR lapply(m@cnms,length) -- but we should have an extractor for this
  nlevel<-sapply(m@flist, function(x) length(levels(x))) #how many inner Blockks in Sigma_i
  ## works as is -- but we should have an extractor
  Vr <- matrix(0,NCOL(Z),NCOL(Z)) #Cov(RanEf)/Var(Error)
  from <- 1
  for(i in 1:k)
  {
    ii<-nlevel[i]
    inner.Blockk<-as.matrix(Sigma.l[[i]])
    to<-from-1+ii*NCOL(inner.Blockk)
    Vr[from:to,from:to]<- inner.Blockk %x% diag(ii)
    from<-to+1
  }
  return(list(
    Vr=Vr, #Cov(RanEf)/Var(Error)
    X=X,
    Z=Z,
    sigmasq=resvar,
    lambda=unique(diag(Vr)),
    y=model.response(model.frame(m)),
    k=k
  ))
}

