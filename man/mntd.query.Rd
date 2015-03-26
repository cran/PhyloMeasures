\name{mntd.query}
\alias{mntd.query}
\title{
Computes the (standardised) value of the Mean Nearest Taxon Distance measure
}
\description{
Calculates the Mean Nearest Taxon Distance (MNTD) for sets of tips on a phylogeny. The same function can also calculate the standardised value of this measure for the given tip sets (the standardised value is equal to minus one times the so-called Nearest Taxon Index, NTI). Note: if the input argument is.standardised is set to TRUE, then this function returns a result only if the input tree is ultrametric.}
\usage{
mntd.query(tree, matrix, is.standardised = FALSE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                must match a tip label on the input tree. If not all values in the matrix are binary, 
                we consider two cases; if the matrix contains only non-negative values, all values are 
                coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                negative value, the function throws an error.}
  \item{is.standardised}{Specifies whether the function returns the standardized value of the MNTD, which reflects the variation in species richness. For each sample R, the MNTD value is standardised by subtracting the mean MNTD and dividing by the standard deviation. The mean and standard deviation are calculated among all tip sets that have the same number of elements as set R, the tip set whose value we want to standardise (default = FALSE)}
}
\value{A vector with the (standardised) MNTD values for each row in the matrix}
\references{Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{mntd.moments}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate mntd values for each community
mntd.query(bird.families,comm)

#Calculate standardised versions
mntd.query(bird.families,comm,TRUE)
}

