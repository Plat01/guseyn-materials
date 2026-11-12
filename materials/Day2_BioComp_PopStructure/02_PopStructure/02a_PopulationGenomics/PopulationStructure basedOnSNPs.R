###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ### 
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###            Day 2: Computing and visualising  population structure from SNP data         ###   
###                                          ----                                           ###
###                    R Script edited by Pr Laurent Gentzbittel, INP-ENSAT                 ###
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
#############################################################################################################


# Loading of the R packages needed for the analysis.
library(adegenet)
library(ape) 
library(parallel)
library(ggplot2)

library(mapplots)  # Data visualisation on maps
library(maptools)  # Tools for reading and handling spatial objects
library(rworldmap) # Mapping global data, vector and raster
library(rworldxtra) #To get the highest resolution
###################
# Data import in R
###################

## -------------------  Using the results of Admixture ------------------- ##

## visualizing cross-valisation analyses : a way to choose among the K's

CVAdmixture <- read.table("CV_sampled.csv", sep = ";", header = TRUE)


x11()
(CV <- ggplot(CVAdmixture, aes(x = K, y = CV)) +
        geom_line( color = "firebrick", size = 1.2) +
        scale_x_continuous(breaks = seq(0,15,1))
    )



## load admixture components for K=5 populations
tblAdmixture <- read.table("KrAllSubset240K-admixture.5.Q")
tblAdmixture$Accession <-  read.table("KrAllSubset240K-admixture.names")[1]
head(tblAdmixture)


ORDER <- order(tblAdmixture[,1], tblAdmixture[,2],tblAdmixture[,3], tblAdmixture[,4], tblAdmixture[,5], decreasing = TRUE)

## syntax is not friendly, because we need both to transpose the data AND order by component.
x11()
barplot(t((tblAdmixture[ORDER, 1:5])), col = funky(5), border = NA,
        xlab = "accessions ", ylab = "Ancestry", 
        names.arg = t(tblAdmixture[ORDER ,6]),
        cex.names = 0.30, las = 2, cex.axis = 1,
        main='Mtr accessions, Admixture analysis'
        )





## -------------------  Using multivariate methods --- ------------------- ##




# data were formated using plink1.9 in plink's raw format
#
## caution :  this (not a very) large dataset that will require lots of computation time to do the course. Use it as your own risks !
## SNPs <- read.PLINK("../data/KrAllSubset240K.raw", map.file = "../data/KrAllSubset240K.map", n.cores = 6, chunkSize = 50)
#
## This i a lighter dataset, for fast computations
 SNPs <- read.PLINK("KrAllSubset16K.raw", map.file = "KrAllSubset16K.map", n.cores = 1)
# chunkSize is the number of genomes to be read at a time.

SNPs


# additional information on accessions 
# shall be in the SAME order as in SNPs
Infos <- read.table("AccessionsGeography.csv", sep=";", header=TRUE)

# geo-cordinates (lat/long) in the 'other' slot of SNPs
other(SNPs)$xy <- Infos[,c("Long", "Lat")]

SNPs
SNPs@other$xy
    

# faster to reload than generating the object again :
save(SNPs, file = "SNPs.RData")
## # if needed reload SNPs object by :
## load("SNPs.RData")



## One can easily derive the distributiong of allele frequencies using:
## glMean: computes the mean of second alleles, i.e. second allele frequencies for each SNP.

myFreq <- glMean(SNPs)

x11()
hist(myFreq, proba = TRUE, col = "gold", xlab = "Allele frequencies",
main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8, lwd = 3)


## However, the distribution of allele frequencies may be more interpretable by restoring its native symmetry:

myFreq <- c(myFreq, 1-myFreq)

x11()
hist(myFreq, proba = TRUE, col = "darkseagreen3", xlab = "Allele frequencies",
main = "Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2, lwd=3)



## ------------------ hierarchical clustering ------------------ ##

D <- dist(as.matrix(SNPs))

## save distance matrix 
save(D, file = "DistanceMatrix.RData")
## reload, without recomputing,  using
## load("DistanceMatrix.RData")


h1 <- hclust(D, method="complete")
h1

x11()
plot(h1, labels = FALSE)
# looks 5 or 6 groups

grp1 <- cutree(h1, k = 5)
head(grp1, 50)


grp2 <- cutree(h1, k = 6)
head(grp2, 100)

## tabulate groups :
table.value(table(grp2, grp1), col.lab=paste("grp", 1:15))
# K = 5?



## ------------------ principal component Analysis ------------------ ##


### PCA , a first, usefull dimension reduction
ptm <- proc.time()

PCA1 <- glPca(SNPs, center = TRUE, scale = FALSE, nf = 200, loadings = TRUE, 
              alleleAsUnit = FALSE, useC = FALSE, parallel = TRUE,
              n.cores = 1, returnDotProd = FALSE, matDotProd = NULL)

(Temps.glPca <- proc.time() - ptm)
## about 5 minutes

## save PCA analysis
save(PCA1, file = "PCA1_glPca.RData")
## reload using
## load("PCA1_glPca.RData")

## save ALL workspace before continuying, in case of problems
save.image(file = "PopStructStep1.RData")


x11()
scatter(PCA1, posi="bottomleft")
title("PCA of Mtr accessions\n axes 1-2, 16.000 SNPs")
## the PCA show


## comparisons with distance-based trees :


tre <- nj(D)
tre
# useful ?
x11()
plot(tre, typ="fan", cex = 0.7)
title("NJ tree of Mtr accessions")


# compare PCA and NJ tree with colors.
# Use mix of RGB of three fist components to create colors. see ?colorplot
x11()
myCol <- colorplot(PCA1$scores, PCA1$scores, transp = TRUE, cex = 4)
abline(h = 0,v = 0, col = "grey")
##add.scatter.eig(PCA1$eig[1:40],2,1,2, posi="bottomright", inset=.05, ratio=.3)
x11()
plot(tre, type = "fan", show.tip = FALSE)
tiplabels(pch = 20, col = myCol, cex = 3.5)
title("NJ tree of Mtr accessions")


## Kmeans :
grp1 <- find.clusters(SNPs, max.n.clust = 40, n.pca = 200,
                      choose.n.clust = FALSE, criterion = "min",
                      n.start = 120, glPca = PCA1)
grp1

# to draw using ggplot
bic <- data.frame(K = seq(1:40), bic = grp1$Kstat)

x11()
(BIC1 <- ggplot(bic, aes(x = K, y = bic)) +
        geom_line( color = "firebrick", size = 1.2))


## DAPC
## We want to assess the relationships between these groups using DAPC.
## Perform the DAPC and store the results in a new object called dapc:

dapc1 <- dapc(SNPs, pop = grp1$grp, scale = FALSE, n.pca = 200, n.da=8,
              glPca = PCA1)


save( dapc1, file= "dapc1_ana.RData")
## load("dapc1_ana.RData")
save.image(file = "PopStructStep2.RData")
## load("dapc1_ana.RData")

#png("dapc1_optim.png")
x11()
scatter(dapc1, 
#        scree.pca = FALSE,
        posi.da="topright",
        pch=20, bg="white", solid=.7,
        col=funky(5)#,
#        col = deepseasun(5)
#        mstree= TRUE
        )
#dev.off()

x11()
scatter(dapc1, scree.da = FALSE, bg="white",
        pch=20, cell=0, cstar=0, col=funky(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))


x11()
scatter(dapc1, xax=2, yax = 3,
        scree.da = FALSE, bg="white",
        pch=20, cell=0, cstar=0, col=funky(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))


# How the dapc is discriminating the groups ?
x11()
par(mfrow=c(3,1))
scatter(dapc1, 1, 1,
        col= funky(5), bg = "white",
        scree.da = FALSE, legend=TRUE, solid=.4)
scatter(dapc1, 2, 2,
        col = funky(5), bg = "white",
        scree.da = FALSE, legend = TRUE, solid=.4)
scatter(dapc1, 3, 3,
        col = funky(5), bg = "white",
        scree.da = FALSE, legend = TRUE, solid=.4)



##  Membership probabilities are based on the retained discriminant functions.
## They are stored in dapc objects in the slot posterior. this may be used to GWAS analyses

class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),3)


## The below figure is the simple graphical translation of the posterior table above. Heat colors
## represent membership probabilities (red=1, white=0); blue crosses represent the prior
## cluster provided to DAPC. Here in most individuals, DAPC classification is consistent

x11()
assignplot(dapc1, subset=1:60)


# Note that this information can also be plotted in a STRUCTURE-like (!) way :
x11()
compoplot(dapc1, posi = "bottomright",
          txt.leg = paste("Cluster", 1:5), lab = "",
          ncol = 1, xlab = "individuals")


# We can also have a closer look at a subset of individuals; for instance, for the first 50 individuals:
x11()
compoplot(dapc1, subset = 1:50, posi="bottomright",
          txt.leg = paste("Cluster", 1:5), lab="",
          ncol = 2, xlab = "individuals")


## drawing a Structure-like plot (need some tricks)
tbl <- dapc1$posterior
tbl <- as.data.frame(tbl)
# major group for each accession
maxi <- function(x){which(x==max(x))}
tbl$Major <- apply( tbl, 1, maxi)
# ordering accessions by group
ORDER <- order(tbl$Major)


## M. truncatula putative structure based on DAPC analysis"

x11()
barplot( t(tbl[ORDER, -6]) , col = funky(5), border = NA,
        xlab = "accessions ", ylab = "Posterior probability", 
        names.arg = rownames(tbl),
        cex.names = 0.30, las = 2, cex.axis = 1,
        main='Mtr accessions, DAPC analysis'
        )



#'
#' # Mapping the first DAPCs onto the geographic space.
#'

# subset SNPs genlight object for accessions that have a geo pos.
which(is.na(other(SNPs)$xy))
which(indNames(SNPs)=="HM288")  # remove "Germany" accession

Aretirer <- c(which(is.na(other(SNPs)$xy)), which(indNames(SNPs)=="HM288") )

## subsetting the SNP dataset :
SNPs2 <- SNPs[ - Aretirer, ]
# checking
which(is.na(other(SNPs2)$xy))


## geographical map :
newmap <- getMap(resolution = "high")

x11(width=15, height=8)
## CairoPDF(file = "./DAPC_BassinMedit.pdf",
##          width = 20, height = 14, onefile = TRUE, family = "Helvetica",
##          title = "R Graphics Output", fonts = NULL, version = "1.1",
##          pointsize=18,paper="A4")
plot(newmap, axes=TRUE, ylim=c(29,45),xlim=c(-10,35), col="gray95")

points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - Aretirer ]])

title("Plot Mtr pop. structure using DAPC", col.main = "firebrick", line = -2, cex = 3)

##dev.off()

