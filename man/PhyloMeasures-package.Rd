\name{PhyloMeasures-package}
\alias{PhyloMeasures-package}
\alias{PhyloMeasures}
\docType{package}
\title{
PhyloMeasures: Fast Computations of Phylogenetic 
Biodiversity Measures
}
\description{ The package provides functions for computing the value and the statistical moments
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
Version: \tab 1.1\cr
Date: \tab 2015-03-20\cr
License: \tab GPL-3 \cr
}
The package contains two types of functions; functions that compute the value of a phylogenetic biodiversity measure, and functions that compute the statistical moments of a measure. We refer to the functions of the first type as query functions, and we refer to the functions of the latter type as moment functions. 

The package provides one query function for each supported phylogenetic measure. The only exception is the CDNT measure, for which three query functions are provided, each handling a different version of this measure. For single-sample measures (these are the PD, the MPD, the MNTD, and the CAC), the corresponding query function receives a tree object (built using package ape), one matrix which represents the tip sets that we want to process, and a logical value that specifies whether for each input sample we want to compute the original value of a measure, or its standardised value. The standardised value of a measure M for a given tip set S is equal to the original value of M for S minus the mean value of M, and this divided by the standard deviation of M for these sets. The mean and the standard deviation values that are used here are computed among all subsets of tips in the input tree that have exactly the same number of elements as S. Query functions of two-sample measures may receive an additional matrix, and also arguments that specify which pairs of input samples should be processed. The query functions are: pd.query, mpd.query, mntd.query, cac.query, cbl.query, cd.query, cdnt.query, cdnt.directed.query, cdnt.averaged.query, phylosor.query, and unifrac.query .

The package also provides one moment function for each supported phylogenetic measure, except the Community Distance Nearest Taxon (CDNT), the Phylogenetic Sorensen's Similarity (PhyloSor), and the Unique Fraction (UniFrac). Each moment function receives a tree object (built using package ape), a vector with non-negative integers that indicates the samples sizes for which we want to compute the statistical moments of the measure, and further arguments that specify the statistical moments that we want to compute for each sample size. Especially for the CAC, the corresponding moment function can compute the k order moment of this measure for any positive integer k. For the rest of the measures, the corresponding moment functions can compute the expectation and the standard deviation for a given sample size. The moment functions provided in the package are: pd.moments, mpd.moments, mntd.moments, cac.moments, cbl.moments, and cd.moments.
}
\author{ Constantinos Tsirogiannis and Brody Sandel

Maintainer: Constantinos Tsirogiannis <analekta@gmail.com>
}
\references{
Faith, D.P. 1992. Conservation evaluation and phylogenetic diversity. Biological Conservation 61: 1-10.

Graham, C.H. and P.V.A. Fine. 2008. Phylogenetic beta diversity: linking ecological and evolutionary processes across space and time. Ecology Letters 11: 1265:1277. 

Lozupone C. and R. Knight. 2005. UniFrac: a New Phylogenetic Method for Comparing Microbial Communities. Applied and Environmental Microbiology,71(12):8228-35.

Swenson, N.G. 2011. Phylogenetic beta diversity metrics, trait evolution and inferring functional beta diversity of communities. PLoS ONE: 6: e21264.

Tsirogiannis, C. and B. Sandel. In prep. Fast computation of measures of phylogenetic beta diversity.

Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient computation of popular phylogenetic tree measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.
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

#Use query function to calculate standardised versions
pd.query(bird.families,comm,TRUE)
}
