\name{pd.moments}
\alias{pd.moments}
\title{
Computes the moments of the Phylogenetic Diversity measure
}
\description{
Calculates the mean and standard deviation of the unrooted Phylogenetic Diversity (PD) for a tree and a vector of tip set sizes.
}
\usage{
pd.moments(tree, sample.sizes, comp.expectation = TRUE, 
           comp.deviation = TRUE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{comp.expectation}{Specifies whether the function should compute the mean (default = TRUE)}
  \item{comp.deviation}{Specifies whether the function should compute the standard deviation (default = TRUE)}
}
\value{If both comp.expectation and comp.deviation are TRUE, the function returns a two-column matrix with one row per element in sample.sizes, where the first column stores the mean PD and the second column stores the standard deviation for this sample size. If only one of comp.expectation or comp.deviation are TRUE, the function returns a vector with the corresponding values instead.}
\references{
Faith, D.P. 1992. Conservation evaluation and phylogenetic diversity. Biological Conservation 61: 1-10.

Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient computation of popular phylogenetic tree measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}

\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{pd.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

pd.moments(bird.families,1:100)
}
