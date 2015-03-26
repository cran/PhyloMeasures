\name{cdnt.directed.query}
\alias{cdnt.directed.query}
\title{Computes the value of the directed Community Distance Nearest Taxon measure}
\description{
Calculates the directed Community Distance Nearest Taxon (dCDNT) given paired sets of tips on a phylogeny. The dCDNT is a beta diversity version of MNTD, giving the mean distance on the tree of tips in a set A to their nearest neighbors in a set B. Note that the measure is asymmetrical, as for taxon t in A the nearest neighbor in B may be taxon s, but for taxon s the nearest neighbor in A might be a third taxon u. Therefore, for every input pair of tip sets A and B, the function returns two values; the dCDNT value from A to B, and the dCDNT value from B to A. The version of this function that computes  the standardised value of the measure is not yet available.
}
\usage{
cdnt.directed.query(tree, matrix.a, matrix.b = NULL, 
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
                      for which the function computes the dCDNT values. Each row in query.matrix
                      indicates a pair of tip sets for which we want to compute the two dCDNT values.
                      Let k and r be the values that are stored in the i-th row of query.matrix, where k
                      is the value stored in the first column and r is the value stored 
                      in the second column. If matrix.b is given, the function computes 
                      both the dCDNT value from the tip set in the k-th row of matrix.a to the tip set
                      in the r-th row of matrix.b, and the dCDNT value from the tip set in the r-th row of matrix.b 
                      to the tip set in the k-th row of matrix.a . If matrix.b is not given,
                      the two values computed are the dCDNT from the k-th to the r-th set of matrix.a,
                      and the dCDNT from the r-th to the k-th set of matrix.a  (default = NULL)}
}
\details{Queries can be given in four ways. If neither matrix.b nor query.matrix are given, the function computes the two dCDNT values for each pair of rows in matrix.a . If matrix.b is given but not query.matrix, the function computes the two dCDNT values for each combination of a row in matrix.a with a row in matrix.b . If query.matrix is given and matrix.b is not, the function returns two dCDNT values for each pair of rows in matrix.a specified by query.matrix. If query.matrix and matrix.b are both given, two dCDNT values are computed for each pair of a row in matrix.a specified by the fist column of query.matrix with a row in matrix.b specified in the second column of query.matrix.}
\value{The dCDNT values for the requested pairs of tip sets. If query.matrix is provided, then the values are returned in a list that contains two elements; each element is a one-dimensional vector storing dCDNT values. Let k and r be the values stored in the i-th row of query.matrix. For the first vector in the returned list, the i-th element of this vector is the dCDNT value from the tip set of the k-th row in matrix.a to the tip set of the r-th row in matrix.a (if matrix.b is not specified), or to the tip set of the r-th row of matrix.b (if this matrix is specified). For the second vector of the output list, the i-th element of this vector is the dCDNT value from the tip set of the r-th row in matrix.a (if matrix.b is not specified), or from the tip set of the r-th row of matrix.b (if this matrix is specified) to the tip set of the k-th row in matrix.a. 

If query.matrix is not provided and matrix.b is provided, the function returns a list which consists of two matrix objects; for the first matrix in this list, entry [i,j] stores the dCDNT value from the tip set specified on the i-th row of matrix.a to the tip set in the j-th row of matrix.b. For the second matrix in the output list, entry [i,j] stores the dCDNT value from the tip set specified on the j-th row of matrix.b to the tip set in the i-th row of matrix.a.

If neither query.matrix nor matrix.b are provided, the function returns a matrix object such that matrix entry [i,j] stores the dCDNT value from the tip set specified on the i-th row of matrix.a to the tip set in the j-th row of this matrix.
}
\references{Graham, C.H. and P.V.A. Fine. 2008. Phylogenetic beta diversity: linking ecological and evolutionary processes across space and time. Ecology Letters 11: 1265:1277. 

Swenson, N.G. 2011. Phylogenetic beta diversity metrics, trait evolution and inferring functional beta diversity of communities. PLoS ONE: 6: e21264.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}
\seealso{\code{\link{cdnt.query}}, \code{\link{cdnt.averaged.query}}}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 10 random communities with 50 families each
comm = matrix(0,nrow = 10,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate all pairwise dCDNT values for communities in comm
cdnt.directed.query(bird.families,comm)
#Note that this matrix is asymmetrical, compare to:
cdnt.query(bird.families,comm)
cdnt.averaged.query(bird.families,comm)

#Calculate pairwise distances from 
#the first two rows of comm to all rows
cdnt.directed.query(bird.families, comm[1:2,],comm)

#Calculate the distances from the first two rows 
#to all rows using the query matrix
qm = expand.grid(1:2,1:10)
cdnt.directed.query(bird.families,comm,query.matrix = qm)   
}
