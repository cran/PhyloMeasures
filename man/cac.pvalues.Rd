\name{cac.pvalues}
\alias{cac.pvalues}
\title{Computes the p-values of the Core Ancestor Cost measure}
\description{
Calculates the p-values of the core ancestor cost (CAC) measure for sets of tips on a phylogeny. The p-values can be calculated under two different null models which maintain species richness. The CAC is calculated by identifying the node in the tree that is the most recent common ancestor of at least chi proportion of the tips in the set, where chi is an input parameter larger than 0.5. The CAC is the distance of this node from the root of the tree.
}
\usage{
cac.pvalues(tree, matrix, chi, null.model="uniform", 
            abundance.weights, reps=1000, seed)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix must match a tip label on the input tree. If not all values in the matrix are binary, we consider two cases; if the matrix contains only non-negative values, all values are coerced to binary ones and a warning message is printed. If the matrix contains at least one negative value, the function throws an error.}
  \item{chi}{A number in the interval (0.5,1]}
  \item{null.model}{A character vector (string) that defines which null model is used for computing the p-values of the measure. There are two possible null models that can be used for computing the p-values: these are "uniform" and "sequential". Both models maintain species richness. More specifically, the available models are defined as follows: 

\itemize{
\item \bold{"uniform"} considers samples
with equal (uniform) probability among all possible tip samples of the same richness.

\item \bold{"sequential"} is an abundance-weighted null model where species samples are chosen based on the same method as R's \emph{sample} function. Unlike the other model (which is computed analytically), this model uses Monte-Carlo randomization. 
}

This argument is optional, and its default value is set to "uniform".
}
  \item{abundance.weights}{A vector of positive numeric values. These are the abundance weights that will be used if option "sequential" is selected. The names stored at the vector must match the names of the tips in the tree. This argument is redundant if the "uniform" model is selected.}
  \item{reps}{ An integer that defines the number of Monte-Carlo random repetitions that will be performed when using the "sequential" model. This argument is redundant if the other null model is selected.}
  \item{seed}{ A positive integer that defines the random seed used in the Monte-Carlo randomizations of the "sequential" model. This argument is optional, and becomes redundant if the other null model is selected.}
}
\value{A vector which stores the computed CAC p-values. The i-th entry in this vector stores the CAC p-value of the i-th row in the input matrix.}
\references{
Tsirogiannis, C. and B. Sandel. 2015. PhyloMeasures: A package for computing phylogenetic
biodiversity measures and their statistical moments. Ecography, doi: 10.1111/ecog.01814, 2015.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{cac.moments}}}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

chi=0.6

#Calculate p-values under the uniform model
cac.pvalues(bird.families,comm,chi, reps=1000)

# Create random abundance weights
weights = runif(length(bird.families$tip.label))
names(weights) = bird.families$tip.label

# Calculate p-values under the sequential model
cac.pvalues(bird.families,comm,chi,null.model="sequential",
            abundance.weights=weights, reps=1000)
}
