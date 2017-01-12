\name{cac.moments}
\alias{cac.moments}
\title{
Computes the statistical moments of the Core Ancestor Cost measure
}
\description{
Calculates the mean and standard deviation of the Core Ancestor Cost (CAC) for a tree and a 
vector of tip set sizes. The means and deviations can be calculated under two different null models which maintain species richness. When the "uniform" null model is used, the function can calculate the first k statistical moments of this measure. The CAC is calculated by identifying the node in the tree that is the most recent common ancestor of at least chi proportion of the tips in the set, where chi is an input parameter larger than 0.5. The CAC is the distance of this node from the root of the tree.
}
\usage{
cac.moments(tree, chi, sample.sizes, k=2,  
            null.model="uniform", abundance.weights, reps=1000, seed)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{chi}{A number in the interval (0.5,1]}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{k}{A positive integer specifying the number of moments to compute (default = 2). If the "sequential" model is selected, the only values that can be used for this argument are either one or two.}

  \item{null.model}{A character vector (string) that defines which null model is used for computing the moments of the measure. There are two possible null models that can be used for computing the moments: these are "uniform" and "sequential". Both models maintain species richness. More specifically, the available models are defined as follows: 

\itemize{
\item \bold{"uniform"} considers samples
with equal (uniform) probability among all possible tip samples of the same richness.

\item \bold{"sequential"} is an abundance-weighted null model where species samples are chosen based on the same method as R's \emph{sample} function. Unlike the other model (which is computed analytically), this model uses Monte-Carlo randomization. 
}

This argument is optional, and its default value is "uniform".
}
  \item{abundance.weights}{A vector of positive numeric values. These are the abundance weights that will be used if option "sequential" is selected. The names stored at the vector must match the names of the tips in the tree. This argument is redundant if the "uniform" model is selected.}
  \item{reps}{ An integer that defines the number of Monte-Carlo random repetitions that will be performed when using the "sequential" model. This argument is redundant if the "uniform" model is selected.}
  \item{seed}{ A positive integer that defines the random seed used in the Monte-Carlo randomizations of the "sequential" model. This argument is optional, and becomes redundant if the "uniform" model is selected.}

}
\value{A k-column matrix with length(sample.sizes) rows. Entry [i,j]  in the matrix gives the j-th moment for the i-th sample size in sample.sizes . The first moment in each row is the mean, and for j larger than one, the j-th returned moment is the raw statistical moment of order j. Given a phylogenetic tree T, a sample size r and a value chi, the 
raw statistical moment of order j for the CAC is defined as: 
\deqn{E_{R \in Sub(T,r)}[\ensuremath{\mathrm{CAC}}_{\chi}(T,R)^k],}{E_{R in Sub(T,r)}[CAC_{chi}(T,R)^k],}
where \eqn{E_{R \in Sub(T,r)}}{E_{R \in Sub(T,r)}} denotes the expectation of a random variable for all tip sets
in T that consist of r tips each.}
\references{Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}

\seealso{\code{\link{cac.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

# Calculate first four raw moments under the uniform model
cac.moments(bird.families,0.75,1:100,k=4)

# Create random abundance weights
weights = runif(length(bird.families$tip.label))
names(weights) = bird.families$tip.label

# Calculate mean and variance under the sequential model
cac.moments(bird.families,0.75,1:100,k=2,
            null.model="sequential", abundance.weights=weights, reps=1000)
}
