\name{cd.moments}
\alias{cd.moments}
\title{
Computes the moments of the Community Distance measure
}
\description{
Calculates the mean and standard deviation of the Community Distance (CD) on a given tree for several pairs of tip set sizes.
}
\usage{
cd.moments(tree, sample.sizes, comp.expectation = TRUE, 
           comp.deviation = TRUE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{sample.sizes}{A two-column matrix of non-negative integers indicating the tip set sizes for which the moments should be calculated. The moments are calculated for each row of the matrix. Let x and y be the values that are stored in the i-th row of sample.sizes. For this row, the CD moments are calculated considering all pairs of tip sets that can be extracted from the input tree such that one set has x elements and the other set has y elements.}
  \item{comp.expectation}{Specifies whether the function returns the mean of the CD (default = TRUE)}
  \item{comp.deviation}{Specifies whether the function returns the standard deviation of the CD (default = TRUE)}

}
\value{If both comp.expectation and comp.deviation are TRUE, the function returns a two-column matrix with one row per element in sample.sizes, such that the first column stores the mean CD and the second column stores the standard deviation of this measure. If only one of comp.expectation or comp.deviation is TRUE, the function returns a vector with the corresponding values instead.}
\references{
Graham, C.H. and P.V.A. Fine. 2008. Phylogenetic beta diversity: linking ecological and evolutionary processes across space and time. Ecology Letters 11: 1265:1277. 

Swenson, N.G. 2011. Phylogenetic beta diversity metrics, trait evolution and inferring functional beta diversity of communities. PLoS ONE: 6: e21264.

Tsirogiannis, C. and B. Sandel. In prep. Fast computation of measures of phylogenetic beta diversity.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{cd.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

cd.moments(bird.families,sample.sizes = expand.grid(1:10,1:10))
}
