\name{mpd.moments}
\alias{mpd.moments}
\title{
Computes the moments of the Mean Pairwise Distance measure
}
\description{
Calculates the mean and standard deviation of the Mean Pairwise Distance (MPD) for a tree and a vector of tip set sizes. The means and deviations can be calculated under three different null models which maintain species richness. 
}
\usage{
mpd.moments(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE, 
            null.model="uniform", abundance.weights, reps=1000, seed)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{comp.expectation}{Specifies whether the function should compute the mean (default = TRUE)}
  \item{comp.deviation}{Specifies whether the function should compute the standard deviation (default = TRUE)}
  \item{null.model}{A character vector (string) that defines which null model is used for computing the moments of the measure. There are three possible null models that can be used for computing the moments: these are "uniform", "frequency.by.richness", and "sequential". All these models maintain species richness. More specifically, the available models are defined as follows: 

\itemize{
\item \bold{"uniform"} considers samples
with equal (uniform) probability among all possible tip samples of the same richness.

\item \bold{"frequency.by.richness"} is an abundance-weighted model where species samples are chosen in a manner similar to the following process; first, each species is selected independently with probability proportional to its abundance. If the resulting sample consists of exactly the same number of elements as the input assemblage then it is used by the null model, otherwise it is tossed and the whole process is repeated.

\item \bold{"sequential"} is an abundance-weighted null model where species samples are chosen based on the same method as R's \emph{sample} function. Unlike the other two models (which are computed analytically), this model uses Monte-Carlo randomization. 
}

This argument is optional, and its default value is set to "uniform".

}
  \item{abundance.weights}{A vector of positive numeric values. These are the abundance weights that will be used if either of the options "frequency.by.richness" or "sequential" are selected. The names stored at the vector must match the names of the tips in the tree. This argument is redundant if the "uniform" model is selected.}
  \item{reps}{ An integer that defines the number of Monte-Carlo random repetitions that will be performed when using the "sequential" model. This argument is redundant if any of the other two null models is selected.}
  \item{seed}{ A positive integer that defines the random seed used in the Monte-Carlo randomizations of the "sequential" model. This argument is optional, and becomes redundant if any of the other two null models is selected.}
}
\value{If both comp.expectation and comp.deviation are TRUE, the function returns a two-column matrix with one row per element in sample.sizes, where the first column stores the mean MPD and the second column stores the standard deviation for this sample size. If only one of comp.expectation or comp.deviation are TRUE, the function returns a vector with the corresponding values instead.}
\references{

Tsirogiannis, C. and B. Sandel. 2015. PhyloMeasures: A package for computing phylogenetic
biodiversity measures and their statistical moments. Ecography, doi: 10.1111/ecog.01814, 2015.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}

\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}

\seealso{\code{\link{mpd.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

# Calculate mean and variance under the uniform model
mpd.moments(bird.families,1:100)

# Create random abundance weights
weights = runif(length(bird.families$tip.label))
names(weights) = bird.families$tip.label

# Calculate mean and variance under the sequential model
mpd.moments(bird.families,1:100,
            null.model="sequential", abundance.weights=weights, reps=1000)
}
