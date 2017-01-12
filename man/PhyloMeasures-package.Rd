\name{PhyloMeasures-package}
\alias{PhyloMeasures-package}
\alias{PhyloMeasures}
\docType{package}
\title{
PhyloMeasures: Fast Computations of Phylogenetic 
Biodiversity Measures
}
\description{ The package supports fast biodiversity computations on large phylogenetic data.
More specifically, the package provides functions for computing the stdandardized values and p-values
of several phylogenetic biodiversity measures. The measures which are supported in the current
version of the package are: the Phylogenetic Diversity (PD), the Mean Pairwise Distance (MPD), 
the Mean Nearest Taxon Distance (MNTD), the Core Ancestor Cost (CAC), the Common Branch Length (CBL), 
the Community Distance (CD), the Community Distance Nearest Taxon (CDNT), the Phylogenetic Sorensen's Similarity (PhyloSor),
and the Unique Fraction (UniFrac).
 }
\details{
\tabular{ll}{
Package: \tab PhyloMeasures\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-12-13\cr
License: \tab GPL-3 \cr
}
The package contains three types of functions; functions that compute the (standardized) value of a phylogenetic biodiversity measure (query functions), functions that compute the p-values of a measure (p-value functions), and functions that compute the statistical moments of a measure (moment functions).  

For the single sample measures, the computation of the standardized values, the statistical moments,
and of the p-values can be done under three different null models. 
We refer to these models as uniform, frequency by richness, and sequential.
All these models maintain species richness, and the latter two models are abundance-weighted. 
The strong advantage of this package is the very efficient algorithms that support the provided functions.
These algorithms make it possible to process very large phylogenies in reasonable time even on a standard
computer. For more details and experimental results, the reader can refer to the main paper for this package,
but also the other works that introduce the efficient algorithms used in the package functions. 
}
\author{ Constantinos Tsirogiannis and Brody Sandel

Maintainer: Constantinos Tsirogiannis <tsirogiannis.c@gmail.com>
}
\references{

Tsirogiannis, C. and B. Sandel. 2015. PhyloMeasures: A Package for Computing Phylogenetic
Biodiversity Measures and Their Statistical Moments. Ecography, doi: 10.1111/ecog.01814, 2015.

Tsirogiannis, C. and B. Sandel. 2016. Fast Computation of Measures of Phylogenetic Beta
Diversity. PLoS ONE, 11(4):  e0151167, doi:10.1371/journal.pone.0151167.

Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient Computation of Popular Phylogenetic Tree Measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New Algorithms for Computing Phylogenetic Biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.

Webb, C.O. 2000. Exploring the Phylogenetic Structure of Ecological Communities: An Example for Rain Forest Trees. The American Naturalist 156: 145-155.
}
\keyword{phylogeny}
\keyword{biodiversity}
\keyword{distance}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Use query function to calculate pd values for each community
pd.query(bird.families,comm)

#Use query function to calculate standardised versions under the uniform model
pd.query(bird.families,comm,TRUE)

# Create random abundance weights
weights = runif(length(bird.families$tip.label))
names(weights) = bird.families$tip.label

#Use query function to calculate standardized versions under the sequential model
pd.query(bird.families,comm,TRUE,null.model="sequential",
         abundance.weights=weights, reps=1000)

}
