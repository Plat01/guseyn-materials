###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                         Day 1: Complete Randomized Factorial Design                     ###
###                                                                                         ###
###         Case Study #1: Variety and Shade effect on lemon tolerance to drought           ### 
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, INP-ENSAT          ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################
#
#The objective is to break down a fixed linear model and to familiarize with R software.

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
# use of the conPloidyt 'sum-to-zero' for ANOVAs
options(contrasts=c('contr.sum','contr.poly'))
#contr.sum: sum of effects of factors = 0

#############################################################################################################

Lemons <- read.table("01_Exmpl1_Citrus.csv", sep = ";", header = TRUE)
str(Lemons)

levels(Lemons$Variety)
levels(Lemons$Shade)


# To impose "Strongshade" as the reference level, it must be in the last position
# By default, in alphabetical order ie "Sun"
Lemons$Shade <- factor(Lemons$Shade, levels = c("Sun", "Shade", "StrongShade"))
levels(Lemons$Shade)


modele <- aov( Index ~ Variety + Shade, data = Lemons)

model.tables(modele, type = "effect")
model.matrix(modele)

summary(modele)

summary(modele, split = list( Shade = list(SunvsStrongShade = 1, ShadevsStrongShade = 2)))

                           
                           
#### understanding SSquares when using matrix notation

## let's fit three successive models :
R0 <- aov( Index ~ 1, data = Lemons)  ## ie only the grand mean
R1 <- aov( Index ~ Variety,  data = Lemons)  ## Grand mean is implicit
R2 <- aov( Index ~ Variety + Shade,  data = Lemons)


## for demonstration, getting residual sum of square of each model :
## residuals : e = Y - X.B
## residuals SSquares = e'.e

(RSS0 <- unlist(summary(R0)[1])[2])   ## cryptic for now. you will learn.
(RSS1 <- unlist(summary(R1)[1])[4])
(RSS2 <- unlist(summary(R2)[1])[6])


## Reduction in Sum of Squares (of residuals) due to factor Variety :
RSS0 - RSS1
## Reduction in Sum of Squares (of residuals) due to factor Shade, after fitting Variety :
RSS1 - RSS2

