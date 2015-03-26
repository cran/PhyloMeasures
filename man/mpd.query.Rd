\name{mpd.query}
\alias{mpd.query}
\title{Computes the (standardised) value of the Mean Pairwise Distance measure}
\description{
Calculates the value of the Mean Pairwise Distance (MPD) measure for sets of tips on a phylogeny. The same function can also calculate the standardised value of this measure for the given tip sets (this is equal to minus one times the so-called Net Relatedness Index, NRI).
}
\usage{
mpd.query(tree, matrix, is.standardised = FALSE)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                must match a tip label on the input tree. If not all values in the matrix are binary, 
                we consider two cases; if the matrix contains only non-negative values, all values are 
                coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                negative value, the function throws an error.}
  \item{is.standardised}{Specifies whether the function should standardize the MPD for variation in species richness. For each tip set R, the observed MPD is standardised by subtracting the mean of the MPD and dividing by the standard deviation of this measure. The mean and standard deviation are calculated among all tip sets that have the same number of elements as set R, the tip set whose value we want to standardise (default = FALSE)}
}
\value{A vector which stores the computed (standardised) MPD values. The i-th entry in this vector stores the MPD value of the i-th row in the input matrix.}
\references{Tsirogiannis, C., B. Sandel and D. Cheliotis. 2012. Efficient computation of popular phylogenetic tree measures. Algorithms in Bioinformatics, LNCS 7534: 30-43.

Webb, C.O. 2000. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. The American Naturalist 156: 145-155.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}
\seealso{\code{\link{mpd.moments}}}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate mpd values for each community
mpd.query(bird.families,comm)

#Calculate standardised versions
mpd.query(bird.families,comm,TRUE)
}

