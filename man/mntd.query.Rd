\name{mntd.query}
\alias{mntd.query}
\title{Computes the (standardized) value of the Mean Nearest Taxon Distance measure}
\description{
Calculates the Mean Nearest Taxon Distance (MNTD) for sets of tips on a phylogeny. The same function can also calculate the standardised value of this measure under three different null models which maintain species richness (the standardised value is equal to minus one times the so-called Nearest Taxon Index, NTI). Note: if the input argument is.standardised is set to TRUE and the null model used is either "uniform" or "frequency.by.richness", then this function returns a result only if the input tree is ultrametric.
}
\usage{
mntd.query(tree, matrix, standardize = FALSE, 
           null.model="uniform", abundance.weights, reps=1000, seed)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix must match a tip label on the input tree. If not all values in the matrix are binary, we consider two cases; if the matrix contains only non-negative values, all values are coerced to binary ones and a warning message is printed. If the matrix contains at least one negative value, the function throws an error.}
  \item{standardize}{Specifies whether the function should standardize the MNTD for variation in species richness. For each tip set S, the observed MNTD is standardized by subtracting the mean MNTD and dividing by the standard deviation of this measure. The mean and standard deviation are calculated among all tip sets that have the same number of elements as set S, the tip set whose value we want to standardize (default = FALSE).}
  \item{null.model}{A character vector (string) that defines which null model is used for computing the standardized values of the measure. There are three possible null models that can be used for computing the standardized values: these are "uniform", "frequency.by.richness", and "sequential". All these models maintain species richness. More specifically, the available models are defined as follows: 

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
\value{A vector which stores the computed (standardized) MNTD values. The i-th entry in this vector stores the (standardized) MNTD value of the i-th row in the input matrix.}
\references{

Tsirogiannis, C. and B. Sandel. 2015. PhyloMeasures: A package for computing phylogenetic
biodiversity measures and their statistical moments. Ecography, doi: 10.1111/ecog.01814, 2015.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{mntd.moments}}}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate mntd values for each community
mntd.query(bird.families,comm)

#Calculate standardized versions under the uniform model
mntd.query(bird.families,comm,TRUE)

# Create random abundance weights
weights = runif(length(bird.families$tip.label))
names(weights) = bird.families$tip.label

#Use query function to calculate standardized versions under the sequential model
mntd.query(bird.families,comm,TRUE,null.model="sequential",
           abundance.weights=weights, reps=1000)
}

