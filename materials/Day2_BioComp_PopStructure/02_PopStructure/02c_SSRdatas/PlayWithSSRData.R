## the tata.xlxs and toto.xlxs are typical excel-like files that contains information about the
## genotypes of a couple of indvidiauls for some SSRs. Please study the format.

## The first goal is to save the tata.xlxs and toto.xlxs files as
## SSRExample1.csv and SSRExample2.csv  files
## using a text (CSV) format.


## reading introductive tutorial for pegas and adegenet packages is higly recommended (google is your best friend)


## ------------------------  loci object ------------------------ ##

library(pegas)

## pegas has the class "loci" which is a simple data frame with a mandatory
## attribute named "locicol" which identifies the columns that are loci; the other
## columns are additional (individual) variables that may be of any kind. The loci
## are coded with factors where the different levels are the observed genotypes with
## the alleles separated with a forward slash, for instance, ‘A/A’ for a classical
## genotype, or ‘132/148’ for a microsatellite locus.
## NO PCR amplication : code 0


## Reading Text Tabular Files
## It is intuitive to organise allelic data in a tabular form where each row is an
## individual and the columns are loci and possibly other variables such as popu-
## lation, spatial locations, and so on.



z <- read.loci("./SSRexample1.csv", loci.sep = ";", col.loci = 2:3, col.pop = 4, row.names = 1)


print(z, details = TRUE)

str(z)

View(z) ## in case of

getAlleles(z)

## We may check that the attribute "locicol" has been set correctly, but usually the user does not need:
attr(z, "locicol")


## Let us have a look on the different options of read.loci:
args(read.loci)



## a more elaborated dataset:
data(jaguar)

print(jaguar, details = TRUE)

str(jaguar)

Fst(jaguar)   ## it also exists a specialized package : ‘mmod’



## ------------------------  genind object ------------------------ ##

library(adegenet)


## adegenet has the class "genind" where individuals are rows and alleles are
## columns in a matrix named tab. This is an S4 class, so the elements are accessed
## with the @ operator (e.g., x@tab). Additional information are stored in other
## slots (@ind.names, @pop, . . . ) The details can be found with class?genind.


toto <- read.table("SSRexample2.csv", header = TRUE, sep = ";", row.names = 1)  ## this indicates that sample names are in the first column
toto

w <- df2genind(toto, ploidy = 2, sep = "/")
w

hw.test(w)

## SNPs data
## In adegenet, SNP data can be handled in two different ways. For relatively small datasets
## (up to a few thousand SNPs) SNPs can be handled as usual codominant markers such as
## microsatellites using genind objects.
## In the case of genome-wide SNP data (from hundreds
## of thousands to millions of SNPs), genind objects are no longer efficient representation
## of the data. In this case, we use genlight objects to store and handle information with
## maximum efficiency and minimum memory requirements.




## converting from ‘loci’ to ‘genind’, package:pegas 
loci2genind(jaguar)


## converting from ‘genind’ to ‘loci’, package:pegas 
as.loci(nancycats)
