\name{cac.query}
\alias{cac.query}
\title{
Computes the (standardised) value of the Core Ancestor Cost measure
}
\description{
Calculates the Core Ancestor Cost (CAC) for sets of tips on a phylogeny. The CAC is calculated by identifying the node in the tree that is the most recent common ancestor of at least chi proportion of the tips in the set. The CAC is the distance of this node from the root of the tree. The same function can also calculate the standardised value of this measure for the given tip sets. 
}
\usage{
cac.query(tree, matrix, chi, is.standardised = FALSE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                  must match a tip label on the input tree. If not all values in the matrix are binary, 
                  we consider two cases; if the matrix contains only non-negative values, all values are 
                  coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                  negative value, the function throws an error.}
  \item{chi}{A number in the interval (0.5,1]}
  \item{is.standardised}{Specifies whether the function should standardise the CAC value according to the variation in species richness. For each input tip set R, the CAC value for R is standardised by subtracting the mean CAC and dividing by the standard deviation of this measure. The mean and standard deviation are calculated among all tip sets that have the same number of elements as set R (default = FALSE)}
}
\value{A vector which stores the computed (standardised) CAC values. The i-th entry in this vector stores the CAC value of the i-th row in the input matrix.}
\references{Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{cac.moments}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate cac values for each community
cac.query(bird.families,comm,0.8)

#Calculate standardised versions
cac.query(bird.families,comm,0.8,TRUE)
}

