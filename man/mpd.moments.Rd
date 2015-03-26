\name{mpd.moments}
\alias{mpd.moments}
\title{
Computes the moments of the Mean Pairwise Distance measure
}
\description{
Calculates the mean and standard deviation of the Mean Pairwise Distance (MPD) for a tree and a vector of tip set sizes.
}
\usage{
mpd.moments(tree, sample.sizes, comp.expectation = TRUE, 
            comp.deviation = TRUE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{comp.expectation}{Specifies whether the function should compute the mean (default = TRUE)}
  \item{comp.deviation}{Specifies whether the function should compute the standard deviation (default = TRUE)}
}
\value{If both comp.expectation and comp.deviation are TRUE, the function returns a two-column matrix with one row per element in sample.sizes, where the first column stores the mean MPD and the second column stores the standard deviation of the MPD for this sample size. If only one of comp.expectation or comp.deviation are TRUE, the function returns a vector with the corresponding values instead.}
\references{Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient computation of popular phylogenetic tree measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{mpd.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

mpd.moments(bird.families,1:100)
}
