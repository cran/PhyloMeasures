\name{cdnt.query}
\alias{cdnt.query}
\title{
Computes the value of the maximised Community Distance Nearest Taxon measure
}
\description{
Calculates the maximised Community Distance Nearest Taxon (mCDNT) given paired sets of tips on a phylogeny. The mCDNT is a beta diversity version of MNTD, and provides the distance on the tree between nearest neighbors in sets A and B. It is derived by calculating  the directed CDNT (dCDNT) from set A to set B, and from set B to set A, and taking the maximum of these two values. A version of this function that computes  the standardised value of the measure is not yet available.
}
\usage{
cdnt.query(tree, matrix.a, matrix.b = NULL, 
           query.matrix = NULL)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix.a}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                  must match a tip label on the input tree. If not all values in the matrix are binary, 
                  we consider two cases; if the matrix contains only non-negative values, all values are 
                  coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                  negative value, the function throws an error}
  \item{matrix.b}{Optional, a second matrix with a similar format as matrix.a}
  \item{query.matrix}{Optional, a two-column matrix specifying the pairs of rows (tip sets)
                      for which the function computes the mCDNT values. Each row in query.matrix
                      indicates a pair of tip sets for which we want to compute the mCDNT value.
                      Let k and r be the values that are stored in the i-th row of query.matrix, where k
                      is the value stored in the first column and r is the value stored 
                      in the second column. If matrix.b is given, the mCDNT value is computed between the 
                      k-th row of matrix.a and the r-th row of matrix.b. If matrix.b is not given,
                      the mCDNT value is computed between the k-th and r-th row of matrix.a (default = NULL)}
}
\details{Queries can be given in four ways. If neither matrix.b nor query.matrix are given, the function computes the mCDNT values for all pairs of rows (tip sets) in matrix.a . If matrix.b is given but not query.matrix, the function computes the mCDNT values for all combinations of a row in matrix.a with rows in matrix.b. If query.matrix is given and matrix.b is not, the function returns the mCDNT values for the pairs of rows in matrix.a specified by query.matrix. If query.matrix and matrix.b are both given, mCDNT values are computed for the rows in matrix.a specified by the first column of query.matrix against the rows in matrix.b specified in the second column of query.matrix}
\value{The mCDNT values for the requested pairs of tip sets. If query.matrix is provided, then the values are returned in an one-dimensional vector. The i-th element of this vector is the mCDNT value for the pair of tip sets 
indicated in the i-th row of query.matrix. If query.matrix is not provided, the mCDNT values are returned in a matrix object;
entry [i,j] in the output matrix stores the mCDNT value between the tip sets specified on the i-th and j-th row of matrix.a (if matrix.b is not specified), or the mCDNT value between the i-th row of matrix.a and the j-th row of matrix.b (if matrix.b is specified)}
\references{Graham, C.H. and P.V.A. Fine. 2008. Phylogenetic beta diversity: linking ecological and evolutionary processes across space and time. Ecology Letters 11: 1265:1277. 

Swenson, N.G. 2011. Phylogenetic beta diversity metrics, trait evolution and inferring functional beta diversity of communities. PLoS ONE: 6: e21264.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{cdnt.directed.query}}, \code{\link{cdnt.averaged.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 10 random communities with 50 families each
comm = matrix(0,nrow = 10,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate all pairwise mCDNT values for communities in comm
cdnt.query(bird.families,comm)

#Calculate pairwise distances from 
#the first two rows of comm to all rows
cdnt.query(bird.families, comm[1:2,],comm)

#Calculate the distances from the first two rows 
#to all rows using the query matrix
qm = expand.grid(1:2,1:10)
cdnt.query(bird.families,comm,query.matrix = qm)
}

