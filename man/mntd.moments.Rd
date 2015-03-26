\name{mntd.moments}
\alias{mntd.moments}
\title{
Computes the moments of the Mean Nearest Taxon Distance measure
}
\description{
Calculates the mean and standard deviation of the Mean Nearest Taxon Distance (MNTD) for a tree and a vector of tip set sizes. Note: this function returns a result only if the input tree is ultrametric.
}
\usage{
mntd.moments(tree, sample.sizes, comp.expectation = TRUE, 
             comp.deviation = TRUE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{comp.expectation}{Specifies whether the function returns the mean of the MNTD (default = TRUE)}
  \item{comp.deviation}{Specifies whether the function returns the standard deviation of the MNTD (default = TRUE)}
}
\value{If both comp.expectation and comp.deviation are TRUE, the function returns a two-column matrix with one row per element in sample.sizes, where the first column stores the mean MNTD and the second column stores the standard deviation for this sample size. If only one of comp.expectation or comp.deviation are TRUE, the function returns a vector with the corresponding values instead.}
\references{Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{mntd.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

mntd.moments(bird.families,1:100)
}
