\name{pd.query}
\alias{pd.query}
\title{Computes the (standardised) value of the unrooted Phylogenetic Diversity measure}
\description{
Calculates the unrooted phylogenetic diversity (PD) measure for sets of tips on a phylogeny. The same function can also calculate the standardised value of this measure for the given tip sets (this is equal to the Phylogenetic Diversity Index, PDI).
}
\usage{
pd.query(tree, matrix, is.standardised = FALSE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                must match a tip label on the input tree. If not all values in the matrix are binary, 
                we consider two cases; if the matrix contains only non-negative values, all values are 
                coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                negative value, the function throws an error.}
  \item{is.standardised}{Specifies whether the function should standardize the PD for variation in species richness. For each tip set R, the observed PD is standardised by subtracting the mean PD and dividing by the standard deviation of this measure. The mean and standard deviation are calculated among all tip sets that have the same number of elements as set R, the tip set whose value we want to standardise (default = FALSE)}
}
\value{A vector which stores the computed (standardised) PD values. The i-th entry in this vector stores the (standardised) PD value of the i-th row in the input matrix}
\references{
Faith, D.P. 1992. Conservation evaluation and phylogenetic diversity. Biological Conservation 61: 1-10.

Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient computation of popular phylogenetic tree measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}
\seealso{\code{\link{pd.moments}}}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate pd values for each community
pd.query(bird.families,comm)

#Calculate standardised versions
pd.query(bird.families,comm,TRUE)
}

