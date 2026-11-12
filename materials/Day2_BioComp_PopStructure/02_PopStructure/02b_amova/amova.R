
## woodmouse                 package:ape                  R Documentation

## Cytochrome b Gene Sequences of Woodmice

## Description:

##      This is a set of 15 sequences of the mitochondrial gene cytochrome
##      _b_ of the woodmouse (_Apodemus sylvaticus_) which is a subset of
##      the data analysed by Michaux et al. (2003). The full data set is
##      available through GenBank (accession numbers AJ511877 to
##      AJ511987).

## Usage:

##      data(woodmouse)
     
## Format:

##      An object of class ‘"DNAbin"’.

## Source:

##      Michaux, J. R., Magnanou, E., Paradis, E., Nieberding, C. and
##      Libois, R. (2003) Mitochondrial phylogeography of the Woodmouse
##      (_Apodemus sylvaticus_) in the Western Palearctic region.
##      _Molecular Ecology_, *12*, 685-697.



library(pegas)


#'
#'
#' Nucleotide diversity and AMOVA
#' from the work of Excoffier et al, 1992, Genetics.
#'
#' AMOVA in **pegas** is multi-levels.


data(woodmouse)
d <- dist.dna(woodmouse, "n")  # simply count the number of differences


### ----------------------  For AMOVA we need to know the populations/groups.
### ----------------------  Can based on sampling, but also on other criteria.


pop1 <- gl(5, 3, labels = LETTERS[1:5])  # artifical populations : 5 pops, 3 mouse per pop

pop1


amova(d ~ pop1)


#' the test is done using permutations.
#' if one need to remove the permutation, use 'nperm=0'
amova(d ~ pop1, nperm = 100) 
#Nothing significant but expected because toy fake dataset

#' Playing by creating continent and pops
conti <- factor(c(rep("A",7), rep("B",8)))
pop2 <- factor(c(rep(1,3),rep(2,4),rep(3,4),rep(4,4)))


amova(d ~ pop2)
amova(d ~ conti/pop2, nperm = 100) #hierarchical model continent/country/city/house...
#Nothing significant but expected because toy fake dataset